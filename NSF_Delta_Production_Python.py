import numpy as np
import math
from scipy import spatial
import matplotlib.pyplot as plt

import sys

class Data(object):
    def __init__(self, geometric_data, temporal_data, spatial_data, model_data):
        self.geometry = geometric_data
        self.temporal = temporal_data
        self.spatial = spatial_data
        self.model_data = model_data
    
    def dxdy(self):
        dx = self.geometry["breath"]/(self.spatial["cols"]-1)
        dy = dx/self.spatial["AR"]
        return dx,dy
    
    def clinoform(self):
        xcols = np.linspace(0, self.geometry["breath"], self.spatial["cols"])
        xtoe_init = (self.geometry["sealevel"] - self.geometry["basement"])/self.geometry["foreslope"]

        tstep = 1
        eta = []
        while True:
            xtoe = xtoe_init + self.geometry["toevel"]*(tstep - 1)*self.temporal["delt"]
           
            xshore = xtoe - (self.geometry["sealevel"] - self.geometry["basement"])/self.geometry["foreslope"]
            temp = np.minimum(np.maximum((xshore-xcols)*self.geometry["topslope"] + self.geometry["sealevel"], self.geometry["sealevel"]),
                                np.maximum((xtoe-xcols)*self.geometry["foreslope"]+ self.geometry["basement"], self.geometry["basement"]))
            eta.append(temp)
            if xtoe >= 0.8*self.geometry["breath"]:
                break
            tstep = tstep + 1

        etabot = np.zeros((tstep, self.spatial["cols"]))
        return (eta, etabot, tstep)

    def init_domain(self):
        eta, etabot , totstep = self.clinoform()
        dx, dy = self.dxdy()
        ytop = eta[i]
        ybot = etabot[i,:]
        ncol = np.floor((ytop - ybot)/dy) + 1
        N = sum(np.floor((ytop - ybot)/dy) + 1)
        phinew = np.max(ytop)*np.ones((N,1))

    
    def simulation(self):
        eta, etabot, totstep = self.clinoform()
        dx, dy = self.dxdy()
        ytop = eta[0]
        ybot = etabot[0]

        ncol = np.round((ytop - ybot)/dy).astype(int) + 1

        N = int(sum(np.round((ytop - ybot)/dy) + 1))
        phinew = np.max(ytop)*np.ones((N,1))
        connew = []
        for jj in range(self.spatial["cols"]):
            if ytop[jj] > self.geometry["sealevel"]:
                connew.append(np.zeros((ncol[jj],1)))
            else:
                connew.append(np.ones((ncol[jj],1)))
        connew = np.concatenate(connew,0)
        # print(connew)
        # print(connew.shape)

        # SIMULATION RUN HERE
        # totstep = 1
        for i in range(totstep):
            print(i)
            ytop = eta[i]
            ybot = etabot[i,:]

            ncolold = ncol
            ncol = np.round((ytop - ybot)/dy).astype(int) + 1
            N = int(sum(np.round((ytop - ybot)/dy) + 1))

            Btop = []
            Btopsea = []
            x = []
            y = []

            cols = self.spatial["cols"]
            for jj in range(cols):
                DeltaX = dx*jj*np.ones((ncol[jj],1))
                DeltaY = np.append(np.linspace(ybot[jj], ybot[jj]+(ncol[jj]-2)*dy, ncol[jj]-1), ytop[jj])
                x.append(DeltaX)
                y.append(DeltaY)
                Btop.append(np.concatenate(y,0).shape[0])
                if ytop[jj] < self.geometry["sealevel"]:
                    Btopsea.append(np.concatenate(y,0).shape[0])
            x = np.concatenate(x,0)
            y = np.concatenate(y,0).reshape(-1,1)
            
            Btop = np.asarray(Btop)
            Btopsea = np.asarray(Btopsea)
            Bbot = np.insert(Btop[:cols-1]+1,0,1)
            # Bbot = np.r_[[1], Btop[:cols-1]+1]

            con = []
            phi = []
            last = 0
            for jj in range(cols):
                first = last
                last = first + ncolold[jj]

                if ncolold[jj] == ncol[jj]:
                    phi.append(phinew[first:last])
                    con.append(connew[first:last])
                
                if ncolold[jj] > ncol[jj]: # errosion
                    phi.append(phinew[first:last-1])
                    con.append(connew[first:last-1])
                
                if ncolold[jj] < ncol[jj]: # deposition
                    # print("======================================")
                    temp_phi = phinew[first:last]
                    temp_phi = np.append(temp_phi,phinew[last]).reshape(-1,1)
                    temp_c = connew[first:last]
                    temp_c = np.append(temp_c, connew[last]).reshape(-1,1)
                    mnode = np.concatenate(phi,0).shape[0]
                    yrat = (y[mnode-1] - y[mnode-2])/(y[mnode-1] - y[mnode-3])
                    idx = len(temp_c)
                    temp_phi[idx-2] = temp_phi[idx-1] - yrat*(temp_phi[idx-1] - temp_phi[idx-3])
                    temp_c[idx-2] = temp_c[idx-1] - yrat*(temp_c[idx-1] - temp_c[idx-3])
                    phi.append(temp_phi)
                    con.append(temp_c)

            phi = np.concatenate(phi,0)
            con = np.concatenate(con,0)
            phinew = phi
            connew = con
            t = []
            first = 0
            for jj in range(cols-1):
                last = first + ncol[jj] + ncol[jj+1]
                xcc = x[first:last]
                ycc = y[first:last]
                points = np.concatenate((xcc,ycc), axis = 1)
                tcc = spatial.Delaunay(points, furthest_site = False, incremental= False, qhull_options=None)
                # print(tcc.simplices)
                # print("-----")
                tcc = tcc.simplices + first
                # print(np.min(tcc))
                t.append(tcc)
                first = first + ncol[jj]
            t = np.concatenate(t,0)
           

            Ntri = t.shape[0]
            N = x.shape[0]
            maxsup = 0
            for ii in range(N):
                # print(np.argwhere(t == ii).shape[0])
                maxsup = max(maxsup, np.argwhere(t == ii).shape[0])
                # print(ii,maxsup)

            maxsup = 2*maxsup
            # maxsup = 20
            # print(t)
            # sys.exit()
            Volp = np.zeros((N,1)) # CV volume
            xmid = np.zeros((Ntri,1)) # element mid point
            ymid = np.zeros((Ntri,1))
            Nx = np.zeros((Ntri,3)) # derivative of shape function
            Ny = np.zeros((Ntri,3))
            tsup = np.zeros((N,1)).astype(int) # triangle in support (count how many elemtent for a given node)
            isup = np.zeros((N,1)).astype(int) # support index location
            sup = np.zeros((N,maxsup)).astype(int) # region of support around a node data point

            for itri in range(Ntri):
                k1 = t[itri,0]
                k2 = t[itri,1]
                k3 = t[itri,2]
                # print(k1,k2,k3)

                v = (x[k2]*y[k3] - x[k3]*y[k2]-x[k1]*y[k3]+x[k1]*y[k2]+y[k1]*x[k3]-y[k1]*x[k2])/2

                Volp[k1] = Volp[k1] + v/3
                Volp[k2] = Volp[k2] + v/3
                Volp[k3] = Volp[k3] + v/3

                xmid[itri] = (x[k1]+x[k2]+x[k3])/3
                ymid[itri] = (y[k1]+y[k2]+y[k3])/3

                Nx[itri,0] = (y[k2]-y[k3])/(2*v)
                Nx[itri,1] = (y[k3]-y[k1])/(2*v)
                Nx[itri,2] = (y[k1]-y[k2])/(2*v)
                Ny[itri,0] = -(x[k2]-x[k3])/(2*v)
                Ny[itri,1] = -(x[k3]-x[k1])/(2*v)
                Ny[itri,2] = -(x[k1]-x[k2])/(2*v)

                sup[k1,isup[k1]] = k2
                sup[k1,isup[k1]+1] = k3
                isup[k1] = isup[k1]+2
                tsup[k1] = tsup[k1]+1

                sup[k2,isup[k2]] = k3
                sup[k2,isup[k2]+1] = k1
                isup[k2] = isup[k2]+2
                tsup[k2] = tsup[k2]+1

                sup[k3,isup[k3]] = k1
                sup[k3,isup[k3]+1] = k2
                isup[k3] = isup[k3]+2
                tsup[k3] = tsup[k3]+1
            
            # print(np.max(t, axis = 0))
            # print(isup)
            # print(sup[0:11,])
            # print(t[0:20])
            # plt.triplot(x.reshape(-1),y.reshape(-1),t)
            # plt.show()
            # sys.exit()

            Big = 1e18
            BCh = np.zeros((N,1))
            BCs = np.zeros((N,1))
            BBh = np.zeros((N,1)) # fixed head values
            BBs = np.zeros((N,1)) # fixed solute concentration values

            BCh[Btop-1] = Big
            BBh[Btop-1] = Big*y[Btop-1]

            BBh[Btopsea-1] = Big*(y[Btopsea-1]-(y[Btopsea-1]-self.geometry["sealevel"])*self.geometry["rho_rel"])
            BCs[Btop-1] = Big
            BBs[Btopsea-1] = Big*self.geometry["c_sea"]

            ap = np.zeros((N,1))
            asup = np.zeros((N,maxsup))
            isup = np.zeros((N,1)).astype(int) # support index so it need to be integer!

            BBvar = np.zeros((N,1))

            kx = np.zeros((Ntri,1))
            ky = np.zeros((Ntri,1))

            for itri in range(Ntri):
                kx[itri] = self.model_data["kx"]
                ky[itri] = self.model_data["ky"]

                cyc = [[0,1,2],[1,2,0],[2,0,1]]
                for node in range(3):
                    ii = cyc[node][0]
                    jj = cyc[node][1]
                    kk = cyc[node][2]

                    k1 = t[itri,ii]
                    k2 = t[itri,jj]
                    k3 = t[itri,kk]

                    Nx1 = Nx[itri,ii]
                    Nx2 = Nx[itri,jj]
                    Nx3 = Nx[itri,kk]
                    Ny1 = Ny[itri,ii]
                    Ny2 = Ny[itri,jj]
                    Ny3 = Ny[itri,kk]

                    delx = (x[k1] + x[k2] + x[k3])/3 - (x[k1] + x[k2])/2
                    dely = (y[k1] + y[k2] + y[k3])/3 - (y[k1] + y[k2])/2

                    face1_k1 = kx[itri]*Nx1*dely-ky[itri]*Ny1*delx
                    face1_k2 = kx[itri]*Nx2*dely-ky[itri]*Ny2*delx
                    face1_k3 = kx[itri]*Nx3*dely-ky[itri]*Ny3*delx
                    # print(face1_k1, face1_k2, face1_k3)

                    BBvar[k1] = BBvar[k1] - ky[itri]*((self.geometry["rho_rel"]-1)/12)*(5*con[k1]+5*con[k2]+2*con[k3])*delx

                    delx = -(x[k1] + x[k2] + x[k3])/3 + (x[k1] + x[k3])/2
                    dely = -(y[k1] + y[k2] + y[k3])/3 + (y[k1] + y[k3])/2

                    face2_k1 = kx[itri]*Nx1*dely-ky[itri]*Ny1*delx
                    face2_k2 = kx[itri]*Nx2*dely-ky[itri]*Ny2*delx
                    face2_k3 = kx[itri]*Nx3*dely-ky[itri]*Ny3*delx
                    # print(face1_k1+face2_k1)
                    # print(k1,k2,k3)

                    ap[k1] = ap[k1] - face1_k1 - face2_k1
                    asup[k1,isup[k1]] = face1_k2 +face2_k2
                    asup[k1,isup[k1]+1] = face1_k3 + face2_k3
                    isup[k1] = isup[k1]+2

                    BBvar[k1] = BBvar[k1] - ky[itri]*((self.geometry["rho_rel"]-1)/12)*(5*con[k1]+2*con[k2]+5*con[k3])*delx

            conver = 1
            phipre = phinew
            sto = self.geometry["sto"]
            delt = self.temporal["delt"]
            # print(asup[1:50])
            # print(ap)
            if i == 657:
                print(asup.shape)
                print(x.shape)
                print(phinew.shape)

            # sys.exit()
            ii = 0
            while conver > 1e-7:
                ii = ii +1
                # print(i)
                # a = phinew[sup].reshape(N,maxsup)
                # print(sup[0:11,:])
                # print(phinew[0:11,:])
                # print(a)
                # if i == 2:
                #     break
                RHS = np.sum(np.multiply(asup,phinew[sup].reshape(N,maxsup)),1)
                phinew = (sto*np.multiply(Volp,phi)+delt*(RHS.reshape(-1,1)+BBvar)+BBh)/(sto*Volp+BCh+ap*delt)
                # print(max(phinew))
                conver = np.max(abs(phinew - phipre))
                # print(conver)
                phipre = phinew
            phi = phinew
            # print(phi)

            # Solute concentration
            aT = self.model_data["aT"]
            aL = self.model_data["aL"]
            Dmol = self.model_data["Dmol"]

            qx = np.zeros((Ntri,1))
            qy = np.zeros((Ntri,1))
            ap = np.zeros((N,1))
            asup = np.zeros((N,maxsup))
            isup = np.zeros((N,1)).astype(int) # support index so it need to be integer!
            for itri in range(Ntri):
                cyc = [[0,1,2],[1,2,0],[2,0,1]]
                for node in range(3):
                    ii = cyc[node][0]
                    jj = cyc[node][1]
                    kk = cyc[node][2]

                    k1 = t[itri,ii]
                    k2 = t[itri,jj]
                    k3 = t[itri,kk]

                    Nx1 = Nx[itri,ii]
                    Nx2 = Nx[itri,jj]
                    Nx3 = Nx[itri,kk]
                    Ny1 = Ny[itri,ii]
                    Ny2 = Ny[itri,jj]
                    Ny3 = Ny[itri,kk]

                    if node == 0:
                        qxval = -kx[itri]*(Nx1*phi[k1] + Nx2*phi[k2] + Nx3*phi[k3])
                        qyval = -ky[itri]*(Ny1*phi[k1] + Ny2*phi[k2] + Ny3*phi[k3])
                        qxmid = qxval
                        qymid = qyval -ky[itri]*((self.geometry["rho_rel"]-1)/3)*(con[k1] + con[k2] + con[k3])

                        qx[itri] = qxmid
                        qy[itri] = qymid

                        qx2 = qxmid**2
                        qy2 = qymid**2
                        qabs = math.sqrt(qx2+qy2)

                        Dxx = aT*qabs + (aL-aT)*qx2/qabs+Dmol*self.geometry["poros"]
                        Dyy = aT*qabs + (aL-aT)*qy2/qabs+Dmol*self.geometry["poros"]
                        Dxy = (aL-aT)*qxmid*qymid/qabs
                    
                    #Face 1
                    qxface = qxval
                    qyface = qyval - ky[itri]*((self.geometry["rho_rel"]-1)/12)*(5*con[k1]+5*con[k2]+2*con[k3])

                    delx = (x[k1] + x[k2] + x[k3])/3 - (x[k1] + x[k2])/2
                    dely = (y[k1] + y[k2] + y[k3])/3 - (y[k1] + y[k2])/2

                    face1_k1 = (Dxx*Nx1 + Dxy*Ny1)*dely - (Dyy*Ny1+Dxy*Nx1)*delx
                    face1_k2 = (Dxx*Nx2 + Dxy*Ny2)*dely - (Dyy*Ny2+Dxy*Nx2)*delx
                    face1_k3 = (Dxx*Nx3 + Dxy*Ny3)*dely - (Dyy*Ny3+Dxy*Nx3)*delx
                    
                    # Upwind
                    qout = qxface*dely - qyface*delx
                    if qout >= 0:
                        face1_k1 = face1_k1 - qout
                    else:
                        face1_k2 = face1_k2 - qout

                    # Face 2
                    qxface = qxval
                    qyface = qyval - ky[itri]*(self.geometry["rho_rel"] -1)/12*(5*con[k1] + 2*con[k2] + 5*con[k3])

                    delx = -(x[k1] + x[k2] + x[k3])/3 + (x[k1] + x[k2])/2
                    dely = -(y[k1] + y[k2] + y[k3])/3 + (y[k1] + y[k2])/2

                    face2_k1 = (Dxx*Nx1 + Dxy*Ny1)*dely - (Dyy*Ny1+Dxy*Nx1)*delx
                    face2_k2 = (Dxx*Nx2 + Dxy*Ny2)*dely - (Dyy*Ny2+Dxy*Nx2)*delx
                    face2_k3 = (Dxx*Nx3 + Dxy*Ny3)*dely - (Dyy*Ny3+Dxy*Nx3)*delx
                    
                    # Upwind
                    qout = qxface*dely - qyface*delx
                    if qout >= 0:
                        face2_k1 = face2_k1 - qout
                    else:
                        face2_k3 = face2_k3 - qout

                    ap[k1] = ap[k1] - face1_k1 - face2_k1
                    asup[k1, isup[k1]] = face1_k2 + face2_k2
                    asup[k1, isup[k1] + 1] = face1_k3 + face2_k3
                    isup[k1] = isup[k1] + 2
                
            conver = 1
            conpre = connew
            eps = self.geometry["poros"]
            delt = self.temporal["delt"]

            ii = 0
            while conver > 1e-7:
                ii = ii +1
                # print(i)
                # a = phinew[sup].reshape(N,maxsup)
                # print(sup[0:11,:])
                # print(phinew[0:11,:])
                # print(a)
                # if i == 2:
                #     break
                RHS = np.sum(np.multiply(asup,connew[sup].reshape(N,maxsup)),1)
                connew = (eps*np.multiply(Volp,con)+delt*(RHS.reshape(-1,1))+BBs)/(eps*Volp+BCs+ap*delt)
                # print(max(phinew))
                conver = np.max(abs(connew - conpre))
                conpre = connew
            con = connew
            con[con < 1e-20] = 0
            con[con > 1.0] = 1.0

            # print(con[700:720])
            # print(con[0:20])
            # print(phi[0:20])
            # print(min(phi))
            # print(max(phi))
            # print(min(con),max(con))
            # print(max(con))
            # print(len(con[con > 0.5]))
        return con
            
        
        plt.triplot(x.reshape(-1),y.reshape(-1),t)
        plt.show()

                    
























geometric_data = {
    "breath" : 3e5, # breath of domain in meter
    "c_sea" : 1, #salt concentration in sea
    "rho_rel" : 1.025, # relative density of saturated saline = rho_sat/ rho_water 
    "poros" : 0.35, # aquifer porosity 
    "sto" : 0.001, # aquifer storage
    "sealevel" : 300, # sea level setting
    "basement" : 250, # depth of basement
    "foreslope" : 0.05, # fore slope
    "topslope" : 0.0005, # fluvial topset slope
    "toevel" : 0.04, # toe velocity over the basement m/days
}

spatial_data = {
    "cols" : 101, #number of columns
    "AR" : 80, #aspect ratio
}

temporal_data = {
    "delt" : 100, # time step (days)
}

model_data = {
    "kx" : 10, # hydraulic conductivity in x direction
    "ky" : 0.1, # hydraulic conductivity in y direction kx/100
    "Dmol" : 0.00001, #m2/day diffusion coefficient
    "aL" : 50, # dispersiivity
    "aT" : 5, # dispersivity
}

test = Data(geometric_data, temporal_data, spatial_data, model_data)
print(test.dxdy())
x,y,z = test.clinoform()
con = test.simulation()
np.savetxt("concentration.txt",con)
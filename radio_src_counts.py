import numpy as np
import pylab as pl

def read_model(filename):

    infile = open(filename,'r')

    fd=[];nS_all=[];nS_fsrq=[];nS_bll=[];nS_ss=[]
    while True:
        line = infile.readline()
        if not line: break

        if line[0]!='#':
            items = line.split()
            fd.append(np.exp(float(items[0])))
            nS_all.append(float(items[1]))
            nS_fsrq.append(float(items[2]))
            nS_bll.append(float(items[3]))
            nS_ss.append(float(items[4]))

    fd = np.array(fd)
    nS_all = np.array(nS_all)
    nS_fsrq = np.array(nS_fsrq)
    nS_bll = np.array(nS_bll)
    nS_ss = np.array(nS_ss)


    return fd, nS_all, nS_fsrq, nS_bll, nS_ss

# ----------------------------------------------------------------------

def read_data(filename):

    infile = open(filename,'r')

    fd=[];nS=[];err1=[];err2=[]
    while True:
        line = infile.readline()
        if not line: break

        if line[0]!='#':
            items = line.split()
            fd.append(np.exp(float(items[0])))
            nS.append(float(items[1]))
            err1.append(float(items[2]))
            err2.append(float(items[3]))

    fd = np.array(fd)
    nS = np.array(nS)
    err1 = np.array(err1)
    err2 = np.array(err2)

    return fd, nS, err1, err2


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

if __name__ == "__main__":

    pl.subplot(111)

    filename = 'massardi_2010_model_21cm.txt'
    fd,nS,nS_fsrq,nS_bll,nS_ss = read_model(filename)
    pl.plot(fd,nS)
    pl.plot(fd,nS_fsrq,ls=':')
    pl.plot(fd,nS_bll,ls=':')
    pl.plot(fd,nS_ss,ls=':')

    filename = 'massardi_2010_data_21cm.txt'
    fd,nS,err1,err2 = read_data(filename)
    pl.errorbar(fd,nS,yerr=[err1,err2],ls='',fmt='o',ms=2)

    pl.xlabel("S [Jy]")
    pl.ylabel(r"$S^{2.5} n(S)$ [Jy$^{1.5}$/sr]")
    pl.loglog()
    pl.show()

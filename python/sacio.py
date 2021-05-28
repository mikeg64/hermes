import numpy as np





def read_step_sac_ascii(filename,step):
    file = open(filename,'rb')

    #read 5 sac file header lines
    
    #1 opozmf_mhd22    #name line 
    header=file.readline()                                                              
    #2      0  0.00000E+00  2  6 10
    head1=file.readline()
    head1=head1.strip()
    head1col=head1.split()
    
    #3 252 252 252    #2D has 2 values
    head2=file.readline()
    head2=head2.strip()
    head2col=head2.split()
    
    #4  1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00    #2D
    #4  1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00    #3D
    head3=file.readline()
    head3=head3.strip()
    head3col=head3.split()
    
    #5 x y h m1 m2 e b1 b2 eb rhob bg1 bg2   gamma eta   grav1 grav2                   #2D
    #5 x y h m1 m2 m3 e b1 b2 b3 eb rhob bg1 bg2 bg3   gamma eta   grav1 grav2 grav3   #3D
    head4=file.readline()
    head4=head4.strip()
    head4col=head4.split()

    
    #2      0  0.00000E+00  2  6 10
    nits=int(head1col[0])
    time=float(head1col[1])
    ndim=int(head1col[2])
    nvar=int(head1col[3])
    nfields=int(head1col[4])
    
    
    #3 252 252    
    if ndim==2:
        odim=[0,0]
        odim[0]=int(head2col[0])
        odim[1]=int(head2col[1])
    if ndim==2:       
        ntot=odim[0]*odim[1]    
    if ndim==3:
        odim=[0,0,0]
        odim[0]=int(head2col[0])
        odim[1]=int(head2col[1])
        odim[2]=int(head2col[2])
        ntot=odim[0]*odim[1]*odim[2]
    
    modelinfo=(header,nits, time, ndim, nvar, nfields,odim,head3,head4)
    #extract useful information from header lines

    
    
    #first skip the step-1 lines including the header
    for i in range(1,step):
        for il in range(ntot):
            dummy=file.readline()
        #skip the header lines
        for il in range(5):
            dummy=file.readline()
                                                                

    
    if ndim==2:
        alldat=np.zeros((odim[0]*odim[1],ndim+nfields))
    elif ndim==3:
        alldat=np.zeros((odim[0]*odim[1]*odim[2],ndim+nfields))   
    
    #extract components from each line
    count=0
    for line in file:
        if count < ntot:
            line=line.strip()
            columns=line.split()
            for i in range(ndim+nfields):
                alldat[count,i]=float(columns[i])
            count=count+1
    
    #using fortran ordering
    #original sac is fortran and same ordering has been adopted
    if ndim==3:
        alldat=np.reshape(alldat,(odim[0],odim[1],odim[2],nfields+ndim),order='F')
    elif ndim==2:
        alldat=np.reshape(alldat,(odim[0],odim[1],nfields+ndim),order='F')
        
    
    	
    file.close()
    
    return alldat,modelinfo






#filename='../../../configs/hydro/3D_128_spic_asc.ini'

#called as
#from sacio import *
#alldat,modelinfo=read_sac_ascii('../../../configs/hydro/3D_128_spic_asc.ini')


def read_sac_ascii(filename):
    file = open(filename,'rb')
                                                                
    #read 5 sac file header lines
    
    #1 opozmf_mhd22    #name line 
    header=file.readline()                                                              
    #2      0  0.00000E+00  2  6 10
    head1=file.readline()
    head1=head1.strip()
    head1col=head1.split()
    
    #3 252 252 252    #2D has 2 values
    head2=file.readline()
    head2=head2.strip()
    head2col=head2.split()
    
    #4  1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00    #2D
    #4  1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00    #3D
    head3=file.readline()
    head3=head3.strip()
    head3col=head3.split()
    
    #5 x y h m1 m2 e b1 b2 eb rhob bg1 bg2   gamma eta   grav1 grav2                   #2D
    #5 x y h m1 m2 m3 e b1 b2 b3 eb rhob bg1 bg2 bg3   gamma eta   grav1 grav2 grav3   #3D
    head4=file.readline()
    head4=head4.strip()
    head4col=head4.split()

    
    #2      0  0.00000E+00  2  6 10
    nits=int(head1col[0])
    time=float(head1col[1])
    ndim=int(head1col[2])
    nvar=int(head1col[3])
    nfields=int(head1col[4])
    
    
    #3 252 252    
    if ndim==2:
    	odim=[0,0]
    	odim[0]=int(head2col[0])
    	odim[1]=int(head2col[1])
    ntot=odim[0]*odim[1]
    
    if ndim==3:
        odim=[0,0,0]
        odim[0]=int(head2col[0])
        odim[1]=int(head2col[1])
        odim[2]=int(head2col[2])
        ntot=odim[0]*odim[1]*odim[2]
    
    modelinfo=(header,nits, time, ndim, nvar, nfields,odim,head3,head4)
    #extract useful information from header lines
    #print("here")
    #print(odim)
    alldat=np.zeros((odim[0]*odim[1],ndim+nfields))
    if ndim==2:
        alldat=np.zeros((odim[0]*odim[1],ndim+nfields))
    elif ndim==3:
        alldat=np.zeros((odim[0]*odim[1]*odim[2],ndim+nfields))   
    
    #extract components from each line
    count=0
    for line in file:
        if count < ntot:
            line=line.strip()
            columns=line.split()
            for i in range(ndim+nfields):
                alldat[count,i]=float(columns[i])
            count=count+1
    
    #using fortran ordering
    #original sac is fortran and same ordering has been adopted
    if ndim==3:
        alldat=np.reshape(alldat,(odim[0],odim[1],odim[2],nfields+ndim),order='F')
    elif ndim==2:
        alldat=np.reshape(alldat,(odim[0],odim[1],nfields+ndim),order='F')
        
    file.close()
    
    return alldat,modelinfo



#filename='../../../configs/sac_test_asc.ini'

#called as
#alldat,modelinfo=read_sac_ascii('../../../configs/sac_test_asc.ini')



def write_sac_ascii(filename, alldat, modelinfo):
    
    file = open(filename,'wb')    
    
    #this script assumes data has been read using a routine such as sac-read3-ascii.py
    #the following variables are assumed
    #nits
    #time
    #ndim
    #nvar
    #nfields
    
    #dim[2] or dim[3]
    
    #gamma
    #eta
    #grav1
    #grav2
    #grav3
    
    #all data is contained in an array alldat of shape nfields+ndim,dim[0],dim[1]
    
    
    #write header lines
    
    #header='sac_test_asc'
    header=modelinfo[0]
    #modelinfo=(header,nits, time, ndim, nvar, nfields,dim,head3,head4)
    #dim=[128,128]
    #ndim=2
    #nfields=12
    dim=modelinfo[6]
    ndim=modelinfo[3]
    nfields=modelinfo[5]
    time=modelinfo[2]
    nits=modelinfo[1]
    nvar=modelinfo[4]
    
    
    head1=str(nits)+" "+str(time)+" "+str(ndim)+" "+str(nvar)+" "+str(nfields)
    
    if ndim==2:
        head2=str(dim[0])+" "+str(dim[1])
    elif ndim==3:
        head2=str(dim[0])+" "+str(dim[1])+" "+str(dim[2])    
    
    #warning may need to explicityly write the adiabatic parameter and correct gravitational parameters here
    head3="1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00"
    
    if ndim==2:
        head4="x y h m1 m2 e b1 b2 eb rhob bg1 bg2   gamma eta   grav1 grav2"
    elif ndim==3:
        head4="x y z h m1 m2 m3 e b1 b2 b3 eb rhob bg1 bg2 bg3   gamma eta   grav1 grav2 grav3"
    
    file.write(header)
    file.write(head1+"\n")     
    file.write(head2+"\n")
    file.write(head3+"\n")
    file.write(head4+"\n")     
        
    if ndim==3:    
        for i3 in range(dim[2]):
            for i2 in range(dim[1]):
                for i1 in range(dim[0]):
                    line=""
                    for j in range(ndim+nfields):
                        line=line+str(alldat[i1,i2,i3,j])+" "
                    line=line+"\n"
                    file.write(line)
                    
    if ndim==2:    
        for i2 in range(dim[1]):
            for i1 in range(dim[0]):
                line=""
                for j in range(ndim+nfields):              
                    line=line+" "+str(alldat[i1,i2,j])
                line=line+"\n"
                file.write(line)    
                
def read_sac_bin(filename):
    file = open(filename,'rb')
        
    file.seek(0,2)
    eof = file.tell()
    file.seek(0,0)
    
    
    header = file.read(79)
    
    nits = np.fromfile(file,dtype=np.int32,count=1)
    
    time = np.fromfile(file,dtype=np.float64,count=1)
    ndim=np.fromfile(file,dtype=np.int32,count=1)
    nvar=np.fromfile(file,dtype=np.int32,count=1)
    nfields=np.fromfile(file,dtype=np.int32,count=1)
    
    dim = np.fromfile(file,dtype=np.int32,count=ndim)[:ndim]
    
    varbuf = np.fromfile(file,dtype=float,count=7)[:7]
    
    #if ndim=2
    head4 = file.read(79)
    
    #if ndim=3
    head3=''
    for i in range(7):
        head3=head3+str(varbuf[i])
                  
     #typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b} CEV;
    
    if ndim==3:
        alldat=np.fromfile(file,dtype=float,count=(nfields+ndim)*dim[0]*dim[1]*dim[2])[:(nfields+ndim)*dim[0]*dim[1]*dim[2]]
        #if size(alldat)<(nw+ndim)*ndata[0]*ndata[1]*ndata[2]:
        #    alldat=resize(alldat,(nw+ndim)*ndata[0]*ndata[1]*ndata[2])
        alldat=np.reshape(alldat,(nfields+ndim,dim[0],dim[1],dim[2],),'C')   
        
    file.close()
    modelinfo=(header,nits, time, ndim, nvar, nfields,dim,head3,head4)

    return alldat,modelinfo


def write_sac_bin(filename, alldat, modelinfo):
    
    file = open(filename,'wb')    
    
    #this script assumes data has been read using a routine such as sac-read3-ascii.py
    #the following variables are assumed
    #nits
    #time
    #ndim
    #nvar
    #nfields
    
    #dim[2] or dim[3]
    
    #gamma
    #eta
    #grav1
    #grav2
    #grav3
    
    #all data is contained in an array alldat of shape nfields+ndim,dim[0],dim[1]
    
    
    #write header lines
    
    #header='sac_test_asc'
    header=modelinfo[0]
    #modelinfo=(header,nits, time, ndim, nvar, nfields,dim,head3,head4)
    #dim=[128,128]
    #ndim=2
    #nfields=12
    dim=modelinfo[6]
    ndim=modelinfo[3]
    nfields=modelinfo[5]
    time=modelinfo[2]
    nits=modelinfo[1]
    nvar=modelinfo[4]
    print(dim,ndim,nfields,nvar)
    
    
    head1=str(nits)+" "+str(time)+" "+str(ndim)+" "+str(nvar)+" "+str(nfields)+"\n"
    
    if ndim==2:
        head2=str(dim[0])+" "+str(dim[1])+"\n"
    elif ndim==3:
        head2=str(dim[0])+" "+str(dim[1])+" "+str(dim[2])+"\n"    
    
    #warning may need to explicityly write the adiabatic parameter and correct gravitational parameters here
    head3="1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00"+"\n"
    
    if ndim==2:
        head4=b"x y h m1 m2 e b1 b2 eb rhob bg1 bg2   gamma eta   grav1 grav2"+"\n"
    elif ndim==3:
        head4=b"x y z h m1 m2 m3 e b1 b2 b3 eb rhob bg1 bg2 bg3   gamma eta   grav1 grav2 grav3"+"\n"
    
    file.write(header.encode('utf-8'))
    file.write(head1.encode('utf-8'))     
    file.write(head2.encode('utf-8'))
    file.write(head3.encode('utf-8'))
    file.write(head4.encode('utf-8'))     
        
    if ndim==3:    
        for i3 in range(dim[2]):
            for i2 in range(dim[1]):
                for i1 in range(dim[0]):
                    line=""
                    for j in range(ndim+nfields):
                        line=line+str(alldat[j,i1,i2,i3,])
                    line=line+"\n"
                    file.write(line.encode('utf-8'))
                    
    if ndim==2:    
        for i2 in range(dim[1]):
            for i1 in range(dim[0]):
                line=""
                for j in range(ndim+nfields):              
                    line=line+" "+str(alldat[j,i1,i2,])
                line=line+"\n"
                file.write(line.encode('utf-8'))    
                

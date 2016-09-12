#Makefile for StarSimulation IoSimTestSimple

##### User configurable options #####
include make_inputs
#include make_inputs_mpi
### End User configurable options ###








	
#	../../IoGenericSimulationLib.o\
#	../../../IoSimulation/IoGenericSteerSimulation.o\
#	../../../IoSimulation/stdsoap2.o\
#	../../../IoSimulation/soapC.o\
#	../../../IoSimulation/soapClient.o
	
ALLOBJS	=	$(obj)/$(BASICEXOBJS)
ALLBIN	=	$(Bin)/smaug



#MYLIBS	=	libAgentSimulation2.a

#@# Targets follow ---------------------------------

##all:	IoSimTestSimple
spicule4:
	cp ../solarmodels/iosmaugparams.h.spicule4_nob ../include/iosmaugparams.h
	cp ../solarmodels/init_user.cu.spicule ../src/init_user.cu
	cp ../solarmodels/boundary.cu.spicule ../src/boundary.cu
	cp ../solarmodels/usersource.cu.spicule4_nob ../src/usersource.cu
	cp ../solarmodels/initialisation_user.h.spicule4_b ../include/initialisation_user.h
spicule4_mpi:
	cp ../solarmodels/iosmaugparams.h.spicule4_nob_mpi ../include/iosmaugparams.h
	cp ../solarmodels/init_user.cu.spicule ../src/init_user.cu
	cp ../solarmodels/boundary.cu.spicule ../src/boundary.cu
	cp ../solarmodels/usersource.cu.spicule4_nob ../src/usersource.cu
	cp ../solarmodels/initialisation_user.h.spicule4_b ../include/initialisation_user.h
spicule5_mpi:
	cp ../solarmodels/iosmaugparams.h.spicule5_nob_mpi ../include/iosmaugparams.h
	cp ../solarmodels/init_user.cu.spicule ../src/init_user.cu
	cp ../solarmodels/boundary.cu.spicule ../src/boundary.cu
	cp ../solarmodels/usersource.cu.spicule4_nob ../src/usersource.cu
	cp ../solarmodels/initialisation_user.h.spicule5_b_mpi ../include/initialisation_user.h
fluxtube1:
	cp ../models/iosmaugparams.h.fluxtube1 ../include/iosmaugparams.h
	cp ../models/init_user.cu.fluxtube1 ../src/init_user.cu
	cp ../models/boundary.cu.fluxtube1 ../src/boundary.cu
	cp ../models/usersource.cu.fluxtube1 ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

spiculejet1:
	cp ../models/iosmaugparams.h.spiculejet1 ../include/iosmaugparams.h
	cp ../models/init_user.cu.spiculejet1 ../src/init_user.cu
	cp ../models/boundary.cu.spiculejet1 ../src/boundary.cu
	cp ../models/usersource.cu.spiculejet1 ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

spicule1:
	cp ../models/iosmaugparams.h.spicule1 ../include/iosmaugparams.h
	cp ../models/init_user.cu.spicule1 ../src/init_user.cu
	cp ../models/boundary.cu.spicule1 ../src/boundary.cu
	cp ../models/usersource.cu.spicule1 ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
spicule1_mpi:
	cp ../models/iosmaugparams.h.spicule1_mpi ../include/iosmaugparams.h
	cp ../models/init_user.cu.spicule1_mpi ../src/init_user.cu
	cp ../models/boundary.cu.spicule1_mpi ../src/boundary.cu
	cp ../models/usersource.cu.spicule1_mpi ../src/usersource.cu
	cp ../models/initialisation_user.h.spicule1_mpi ../include/initialisation_user.h
sedov:
	cp ../models/iosmaugparams.h.sedov ../include/iosmaugparams.h
	cp ../models/init_user.cu.sedov ../src/init_user.cu
	cp ../models/boundary.cu.sedov ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
sedov3d:
	cp ../models/iosmaugparams.h.sedov3d ../include/iosmaugparams.h
	cp ../models/init_user.cu.sedov3d ../src/init_user.cu
	cp ../models/boundary.cu.sedov3d ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
	cp ../models/initialisation_user.h.sedov3d ../include/initialisation_user.h
bach:
	cp ../models/iosmaugparams.h.bach ../include/iosmaugparams.h
	cp ../models/init_user.cu.bach ../src/init_user.cu
	cp ../models/boundary.cu.bach ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

ah:
	cp ../models/iosmaugparams.h.ah ../include/iosmaugparams.h
	cp ../models/init_user.cu.ah ../src/init_user.cu
	cp ../models/boundary.cu.ah ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

bw:
	cp ../models/iosmaugparams.h.bw ../include/iosmaugparams.h
	cp ../models/init_user.cu.bw ../src/init_user.cu
	cp ../models/boundary.cu.bw ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
rt:
	cp ../models/iosmaugparams.h.rt ../include/iosmaugparams.h
	cp ../models/init_user.cu.rt ../src/init_user.cu
	cp ../models/boundary.cu.rt ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
kh:
	cp ../models/iosmaugparams.h.kh ../include/iosmaugparams.h
	cp ../models/init_user.cu.kh ../src/init_user.cu
	cp ../models/boundary.cu.kh ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
ot:
	cp ../models/iosmaugparams.h.ot ../include/iosmaugparams.h
	cp ../models/init_user.cu.ot ../src/init_user.cu
	cp ../models/boundary.cu.ot ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
ot_mpi:
	cp ../models/iosmaugparams.h.ot_mpi ../include/iosmaugparams.h
	cp ../models/init_user.cu.ot_mpi ../src/init_user.cu
	cp ../models/boundary.cu.ot_mpi ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

kink1:
	cp ../models/iosmaugparams.h.kink1 ../include/iosmaugparams.h
	cp ../models/init_user.cu.kink1 ../src/init_user.cu
	cp ../models/boundary.cu.kink1 ../src/boundary.cu
	cp ../models/usersource.cu.kink1 ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
kink2d:
	cp ../models/iosmaugparams.h.kink2d ../include/iosmaugparams.h
	cp ../models/init_user.cu.kink2d ../src/init_user.cu
	cp ../models/boundary.cu.kink2d ../src/boundary.cu
	cp ../models/usersource.cu.kink2d ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
alf:
	cp ../models/iosmaugparams.h.alf ../include/iosmaugparams.h
	cp ../models/init_user.cu.alf ../src/init_user.cu
	cp ../models/boundary.cu.alf ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
cleanobjs:
	rm -f $(BASICEXOBJS) ../include/*.cuh

cleanbin:
	rm -f $(ALLBIN)


clean:	cleanobjs cleanbin

cleanall:	cleanobjs cleanbin

all: smaug
	$(CC) -g $(LIBS)  $(incDirs) -o $(Bin)/smaug $(ALLOBJS) $(MYLIBS)

allsac:$(ALLOBJS)
	$(NVCC) -g $(LIBS)  $(incDirs) -o $(Bin)/smaug $(ALLOBJS) $(MYLIBS)
#@# Dependency rules follow -----------------------------


#ioshalloww:$(BASICEXOBJS)
#	$(NVCC) --allow-multiple-defintion -g $(LIBS)  $(incDirs) -o $(Bin)/ioshalloww $(BASICEXOBJS) $(MYLIBS)
smaug:$(BASICEXOBJS)
	rm iosmaug.o
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs)  -o ../obj/iosmaug.o -c iosmaug.cpp
	$(CC) $(LINKFLAGS)  $(LIBS)  $(incDirs) -o $(Bin)/smaug $(ALLOBJS) $(MYLIBS)




#IoTest.o: IoTest.f90
#	$(FF) -o $@ -c $<
#step.o: step.cu
#	$(NVCC) -o $@  --device-emulation -c $<
#step.o: step.cu
#	$(NVCC) -o $@ --ptxas-options=-v -c $<
nushk1.o: nushk1.cu
	sed 's/MODID/nshk/g' < gradops.cu > ../include/gradops_nshk.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifvisc1.o: hyperdifvisc1.cu
	sed 's/MODID/hdv1/g' < gradops.cu > ../include/gradops_hdv1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifvisc1l.o: hyperdifvisc1l.cu
	sed 's/MODID/hdv1l/g' < gradops.cu > ../include/gradops_hdv1l.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifvisc1r.o: hyperdifvisc1r.cu
	sed 's/MODID/hdv1r/g' < gradops.cu > ../include/gradops_hdv1r.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifviscmax.o: hyperdifviscmax.cu
	sed 's/MODID/hdvmax/g' < gradops.cu > ../include/gradops_hdvmax.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifrhosource1.o: hyperdifrhosource1.cu
	sed 's/MODID/hdr1/g' < gradops.cu > ../include/gradops_hdr1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifesource1.o: hyperdifesource1.cu
	sed 's/MODID/hde1/g' < gradops.cu > ../include/gradops_hde1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifbsource1.o: hyperdifbsource1.cu
	sed 's/MODID/hdb1/g' < gradops.cu > ../include/gradops_hdb1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifbsourcene1.o: hyperdifbsourcene1.cu
	sed 's/MODID/hdbne1/g' < gradops.cu > ../include/gradops_hdbne1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifmomsource1.o: hyperdifmomsource1.cu
	sed 's/MODID/hdm1/g' < gradops.cu > ../include/gradops_hdm1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifmomsourcene1.o: hyperdifmomsourcene1.cu
	sed 's/MODID/hdmne1/g' < gradops.cu > ../include/gradops_hdmne1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
centdiff1.o: centdiff1.cu
	sed 's/MODID/cd1/g' < gradops.cu > ../include/gradops_cd1.cuh
	sed 's/MODID/cd1/g' < dervfields.cu > ../include/dervfields_cd1.cuh
	sed 's/MODID/cd1/g' < usersource.cu > ../include/usersource_cd1.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
centdiff2.o: centdiff2.cu
	sed 's/MODID/cd2/g' < gradops.cu > ../include/gradops_cd2.cuh
	sed 's/MODID/cd2/g' < dervfields.cu > ../include/dervfields_cd2.cuh
	sed 's/MODID/cd2/g' < usersource.cu > ../include/usersource_cd2.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
computedervfields.o: computedervfields.cu
	sed 's/MODID/cdf/g' < gradops.cu > ../include/gradops_cdf.cuh
	sed 's/MODID/cdf/g' < dervfields.cu > ../include/dervfields_cdf.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS)   -c $<
boundary.o: boundary.cu
	sed 's/MODID/b/g' < gradops.cu > ../include/gradops_b.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
init.o: init.cu
	sed 's/MODID/i/g' < gradops.cu > ../include/gradops_i.cuh
	sed 's/MODID/i/g' < init_user.cu > ../include/init_user_i.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
mpiutil.o: mpiutil.cu
	sed 's/MODID/mpiu/g' < gradops.cu > ../include/gradops_mpiu.cuh
	sed 's/MODID/mpiu/g' < init_user.cu > ../include/init_user_mpiu.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
update.o: update.cu
	sed 's/MODID/u/g' < gradops.cu > ../include/gradops_u.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
advance.o: advance.cu
	sed 's/MODID/adv/g' < gradops.cu > ../include/gradops_adv.cuh
	sed 's/MODID/adv/g' < dervfields.cu > ../include/dervfields_adv.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
divb.o: divb.cu
	sed 's/MODID/db/g' < gradops.cu > ../include/gradops_db.cuh
	sed 's/MODID/db/g' < dervfields.cu > ../include/dervfields_db.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
.cpp.o:
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)  -o $@ -c $<

.c.o:
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)    -o $@ -c $<






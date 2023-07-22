# 3D parallel geometric multigrid makefile
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : January 5th 2016

# Load the common configuration file
include ../../config.mk

iflags=
lflags=-L.

probs=problem_simple.o manu_ps_vec.o helmholtz.o
l_objs=buffer.o common.o geometry.o mat3.o
objs=$(l_objs)
p_src=$(patsubst %.o,%.cc,$(probs))
src=$(patsubst %.o,%.cc,$(objs))
execs=gs_parallel mg_parallel gmap_test factors iir_test lapack_test mg_converge helm_test

all:
	$(MAKE) -C ../../shared
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(mpicxx) $(iflags) -MM $(p_src) $(src) >Makefile.dep

include Makefile.dep

libmg3d.a: $(l_objs)
	rm -f libmg3d.a
	ar rs libmg3d.a $^

helm_test: helm_test.cc helmholtz.o libmg3d.a
	$(mpicxx) $(cflags) $(iflags) -o $@ $(lflags) $^ $(lp_lflags)

mg_converge: mg_converge.cc manu_ps_vec.o libmg3d.a
	$(mpicxx) $(cflags) $(iflags) -o $@ $(lflags) $^ $(lp_lflags)

gs_parallel: gs_parallel.cc problem_simple.o libmg3d.a
	$(mpicxx) $(cflags) $(iflags) -o $@ $(lflags) $^ $(lp_lflags)

mg_parallel: mg_parallel.cc manu_ps_vec.o libmg3d.a
	$(mpicxx) $(cflags) $(iflags) -o $@ $(lflags) $^ $(lp_lflags)

factors: factors.cc
	$(cxx) $(cflags) $(iflags) -o $@ $<

iir_test: iir_test.cc
	$(cxx) $(cflags) $(iflags) -o $@ $<

gmap_test: gmap_test.cc
	$(cxx) $(cflags) $(iflags) -o $@ $<

lapack_test: lapack_test.cc
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lp_lflags)

%.o: %.cc
	$(mpicxx) -Winline $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs) $(probs) libmg3d.a

.PHONY: clean all executables depend

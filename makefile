include /home/junteng/lib/ARPACK/ARmake.inc

PCXX      = mpixlcxx_r
PCXXFLAGS = 

ppca: ppca.o
	$(PFC) $(PFFLAGS) ppca.o $(PLIBS) -o ppca.x

cov: cov.o
	$(PCXX) $(PCXXFLAGS) cov.o -o cov.x

cov.o: cov.cpp
	$(PCXX) -c -o $@ $<

clean:
	rm *.o
	rm *.x

include /home/junteng/lib/arpack/ARPACK/ARmake.inc

ppca: ppca.o
	$(PFC) $(PFFLAGS) ppca.o $(PLIBS) -o ppca.x

cov: 
	$(PCXX) $(PCXXFLAGS) cov.cpp -o cov.x

clean:
	rm *.o
	rm *.x

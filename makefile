include /home/junteng/lib/ARPACK/ARmake.inc

ppca: ppca.o
	$(PFC) $(PFFLAGS) ppca.o $(PLIBS) -o ppca.x

clean:
	rm *.o
	rm *.x

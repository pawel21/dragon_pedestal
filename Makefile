lib: translatezFits.o
	clang translatezFits.o -shared -o translatezFits.so
translatezFits.o: translatezFits.c
	clang -c -fPIC translatezFits.c -o translatezFits.o
clean:
	rm translatezFits.o translatezFits.so

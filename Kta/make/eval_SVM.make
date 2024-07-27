CFLAGS= -O3 -msse2 -mfpmath=sse

eval_SVM: eval_SVM.o algebre.o biblio.o
	gcc eval_SVM.o algebre.o biblio.o -lm $(CFLAGS) -o eval_SVM

algebre.o: algebre.c
	gcc -c algebre.c -o algebre.o $(CFLAGS)

biblio.o: biblio.c
	gcc -c biblio.c -o biblio.o $(CFLAGS)

eval_SVM.o: eval_SVM.c
	gcc -c eval_SVM.c -o eval_SVM.o $(CFLAGS)

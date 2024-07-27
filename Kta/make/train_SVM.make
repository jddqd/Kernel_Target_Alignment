CFLAGS= -O3 -msse2 -mfpmath=sse

train_SVM: train_SVM.o algebre.o
	gcc train_SVM.o algebre.o -lm $(CFLAGS) -o train_SVM

algebre.o: algebre.c
	gcc -c algebre.c -o algebre.o $(CFLAGS)

train_SVM.o: train_SVM.c
	gcc -c train_SVM.c -o train_SVM.o $(CFLAGS)

CFLAGS= -O3 -msse2 -mfpmath=sse

train_SVM: calcul_theta.c
	gcc calcul_theta.c -o calcul_theta $(CFLAGS) -lm

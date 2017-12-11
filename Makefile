#
#
CC = mpicc
OBJ = main.o sorceterm.o FFT_3D_complex.o iteration.o assign_matrix.o err.o Combine_FFT.o label_corr.o corr_procs.o Ini_label.o
solution: $(OBJ)
	$(CC) -o solution $(OBJ) -lfftw3 -lm
main.o: main.c
	gcc -c -g main.c -lfftw3 -lm
sorceterm.o: sorceterm.c
	gcc -c -g sorceterm.c -lfftw3 -lm
FFT_3D_complex.o: FFT_3D_complex.c
	gcc -c -g FFT_3D_complex.c -lfftw3 -lm
iteration.o: iteration.c
	gcc -c -g iteration.c -lfftw3 -lm
assign_matrix.o: assign_matrix.c
	gcc -c -g assign_matrix.c -lfftw3 -lm
err.o: err.c
	gcc -c -g err.c -lfftw3 -lm
Combine_FFT.o: Combine_FFT.c
	gcc -c -g Combine_FFT.c -lfftw3 -lm
label_corr.o: label_corr.c
	gcc -c -g label_corr.c -lfftw3 -lm
corr_procs.o: corr_procs.c
	gcc -c -g corr_procs.c -lfftw3 -lm
Ini_label.o: Ini_label.c 
	gcc -c -g Ini_label.c -lfftw3 -lm

clean:
	rm -f *.o


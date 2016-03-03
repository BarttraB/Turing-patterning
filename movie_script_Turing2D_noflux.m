

L6=load('L_E4_2D_8_noflux.txt');
L6m=matrix_maker2D_3(L6,263);
L6mov=movie_maker5_nox(L6m,1,300,1);
movie2avi(L6mov, 'Turing2Dmovie2_noflux', 'compression', 'None')

%***to convert avi to mpeg
%ffmpeg -i Turing2Dmovie.avi -r 20 Turing2Dmovie.mpg

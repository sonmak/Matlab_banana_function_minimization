clc
clear
close all

%symboliczny graient i macierz hessego

syms x y a b

f = (1 - x + a) ^ 2 + 100 * (y - b - (x - a) ^ 2) ^ 2

hessian(1,1) = diff(f,x,2)
hessian(1,2) = diff(f,x,y)
hessian(2,1) = diff(f,y,x)
hessian(2,2) = diff(f,y,2)

gradient = gradient(f,[x,y])

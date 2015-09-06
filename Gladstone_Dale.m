function [Dens,Dens_av]=Gladstone_Dale(n2,xc,zc)


% Galdstone-Dale constant for saline-water
% https://books.google.co.il/books?id=DJCKI5qQdiAC&pg=PA119&lpg=PA119&dq=gladstone+dale+constant+water&source=bl&ots=q1XkV2tiYM&sig=cdTxBmgIsFh_k0fJrZ-cZVsKA_M&hl=it&sa=X&sqi=2&ved=0CCkQ6AEwAWoVChMI6uXTte2OxgIVUlnbCh14ZwQE#v=onepage&q=gladstone
% search better for that value
G=0.335;  %[g/mL]
S_out=(n2-1)./G;

Dens.x=xc;
Dens.z=zc;
Dens.f=S_out;

Dens_av=mean((S_out));
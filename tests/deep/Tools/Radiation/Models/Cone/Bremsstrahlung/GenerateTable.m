%% clear variables, define r0
clc, %clear
m = 9.10938356*10^(-31);
e = 1.60217662*10^(-19);
c = 2.99792458 *10^8;
epsilon_0 = 8.85418782*10^(-12);
r0 = e^2/(4*pi*epsilon_0*m*c^2); 
%% create table with random values
clc
%nspecies_tab = sort(round(5*rand(1,10)+ones(1,10)));
format long

%Z_tab = round(7*rand(1, sum(nspecies_tab))+ones(1, sum(nspecies_tab)));
%Z0_tab = round(Z_tab.*rand(1, sum(nspecies_tab)));
%density_tab = 10^(19)*(ones(1, sum(nspecies_tab))+10*round(rand(1, sum(nspecies_tab)), 5));
I_tab = power(nspecies_tab, Z_tab, density_tab, r0);

Final_table1 = [nspecies_tab; I_tab]';
Z_fin = Z_tab';
Z0_fin = Z0_tab';
dens_fin = density_tab';
%Final_table2 = [Z; Z0; density]';


%% TEST functions
'new'
Z = [8,5,3,2];
Z0 = 0;
density = [1.07350000000000e+19,7.80040000000000e+19,8.05950000000000e+19,7.45130000000000e+19];
nspecies = [1, 1, 2];
I = power(nspecies, Z, density, r0)

%% MORE TESTS

Z = [3 2];
density = [8.059500e+19, 7.451300e+19];
nspecies = 2;
I = power_light(nspecies, Z, density, r0)

%% Functions

function pow = power(nspecies, Z, density, r0)

I = zeros(1, length(nspecies));
k = 1;
i = 1;
for j=1:length(nspecies)
    while i <= nspecies(j) + k-1
      i
      I(j) = I(j) + density(i)*eval_4BS(Z(i), r0);
      i = i +1;
    end
    k = i;
end
pow = I;
end


function pow = power_light(nspecies, Z, density, r0)

I = 0;

    for i=1:nspecies
      i
      I = I + density(i)*eval_4BS(Z(i), r0);
     end
pow = I;
end

function val = eval_4BS(Z, r0)
    Z2fakt = 4*Z^2*r0^2/137;
    lnfakt = log(183) - 0.5*log(Z)+1/18;
    val = Z2fakt*lnfakt;
end
%% clear variables, define r0
clc, clear
m = 9.10938356*10^(-31);
e = 1.60217662*10^(-19);
c = 2.99792458 *10^8;
epsilon_0 = 8.85418782*10^(-12);
h =  6.62607015*10^(-34);
r0 = e^2/(4*pi*epsilon_0*m*c^2);
alpha = e^2/(2*epsilon_0*h*c);
%% create table with random values - total emission
clc
nspecies_tab = sort(round(5*rand(1,10)+ones(1,10)));
format long

Z_tab = round(7*rand(1, sum(nspecies_tab))+ones(1, sum(nspecies_tab)));
Z0_tab = round(Z_tab.*rand(1, sum(nspecies_tab)));
density_tab = 10^(19)*(ones(1, sum(nspecies_tab))+10*round(rand(1, sum(nspecies_tab)), 5));
I_tab = power(nspecies_tab, Z_tab, density_tab, r0);

Final_table1 = [nspecies_tab; I_tab]';
Z_fin = Z_tab';
Z0_fin = Z0_tab';
dens_fin = density_tab';
%Final_table2 = [Z; Z0; density]';
%% create table of values, generated from before - spectrum
clc
format long

Z_tab = [7 4 4 4 2 3]; %Only one nspecies value (6), still gives 10 referens-values
Z0_tab = [1 0 4 2 2 2];
density_tab = [1.487400e+19 6.527300e+19 3.748100e+19 3.415000e+19 3.431500e+19 2.541600e+19];
gamma = [1 5 10 15 25 30 40 50 75 100];
Integrated_spectrum = zeros(1, numel(gamma));
wl = linspace(1, 50, 50);

for i=1:length(gamma)
    spectrum = Generate_spectrum(Z_tab, Z0_tab, density_tab, wl, gamma(i), r0, alpha);
    Integrated_spectrum(i) = Integrate_spectrum(spectrum, wl);
end

Final_table_spec = [gamma; Integrated_spectrum]';
wl_fin = wl;
%% Table for 3BN
clc; format long;
Z = [2 7 10 15];
dens = [1.649e19 11.467e19 7.3903e19 9.1003e19];
gamma = [5 13 38 67];
Integrated_spectrum = zeros(1, numel(gamma));
wl = linspace(1, 50, 50);


for i=1:length(gamma)
    spectrum = spec_3BN(Z, dens, wl, gamma(i), r0, alpha);
    Integrated_spectrum(i) = Integrate_spectrum(spectrum, wl);
end

Final_table_3BN = [Z; dens; gamma; Integrated_spectrum]';

%% Functions

function val = Get_4BS(Z, r0)
    Z2fakt = 4*Z^2*r0^2*alpha;
    lnfakt = log(183) - 0.5*log(Z)+1/18;
    val = Z2fakt*lnfakt;
end

function spec = Generate_spectrum(Z, Z0, density, k, gamma_in, r0, alpha)
    Nej = Z-Z0;
    A_bar = (9 * pi * Nej.^2).^(1/3)./(2*alpha*Z);
    n_species = length(Z);
    nr_wl = numel(k); %Number of wavelengths
    prefact = 4*r0^2*alpha;
    

    if length(gamma_in) == 1
        gamma_in = gamma_in*ones(1, nr_wl);
    end
    p_in = sqrt(gamma_in.^2-1);
    gamma_out = gamma_in - k; %Final gamma-factor
    
    dsigma_dp = zeros(1, nr_wl);
    for j=1:n_species
        for i=1:nr_wl
            
            if gamma_out(i) < 1
              continue;
            end
            q0 = p_in(i) - sqrt(gamma_out(i)^2-1) - k(i);
        
            sum_part1 = (1+(gamma_out(i)/gamma_in(i))^2)*Integral_1(Z(j), Z0(j), q0, A_bar(j));
            sum_part2 = -2*gamma_out(i)/(3*gamma_in(i))*Integral_2(Z(j), Z0(j), q0, A_bar(j));
            species_cont = density(j)*(sum_part1+sum_part2);
            
            dsigma_dp(i) = dsigma_dp(i) + prefact/(k(i))*species_cont;
        end
        
    end
    spec = dsigma_dp; %sqrt((gamma_out).^2-1)/(gamma_out)* dsigma_dp;
end

function Int = Integral_1(Z, Z0, q0, a_bar)
    Ne = Z-Z0;
    fun = @(q) (Z - Ne/(1+(q*a_bar)^(3/2)))^2*(q-q0)^2/q^3;
    int_val = integral(fun, q0, 1, 'ArrayValued', true);
    Int = Z^2 + int_val;
end

function Int = Integral_2(Z, Z0, q0, a_bar)
    Ne = Z-Z0;
    fun = @(q) (Z - Ne/(1+(q*a_bar)^(3/2)))^2 * (q^3 + 3*q*q0^2* (1-2*q*q0^2*log(q/q0)) - 4*q0^3)/q^4;
    int_val = integral(fun, q0, 1, 'ArrayValued', true);
    Int = 5/6*Z^2 + int_val;
end

function Int_spec = Integrate_spectrum(spec, wl)
    nwl = numel(wl);
    Int = 0.5*(spec(1)+spec(nwl));
    Int = Int + sum(spec(2:nwl-1));
    Int_spec = Int*(wl(2)-wl(1));
end

function spec= spec_3BN(Z, dens, wl, gamma_in, r0, alpha)
    nr_wl = numel(wl);
    d1sigma = zeros(1,nr_wl);
    for i=1:nr_wl
        k = wl(i);
        if k >= gamma_in-1
            continue;
        end
        
        p0 = sqrt(gamma_in.^2 -1);
        E0 = sqrt(1+p0.^2);
        E = E0-k;
        p = sqrt(E.^2-1);

        PreFactor = r0^2*alpha*p./(p0.*k);

        E0E = E0.*E;
        p0p = p0.*p;

        Eps  = 2*log(E+p);
        Eps0 = 2*log(E0+p0);

        L = 2*log( (E0E + p0p - 1)./k );

        Term1 = 4/3 - 2*E0E.* (p.^2 + p0.^2)./(p0p.^2) ...
            + Eps0.*E./p0.^3 + Eps.*E0./p.^3 - Eps.*Eps0./(p0p);

        Term2 = L.*( 8*E0E./(3*p0p) + k.^2./(p0p.^3) .*( E0E.^2 + p0p.^2 ) );

        Term3 = L.* k./(2*p0p) .*( (E0E + p0.^2).*Eps0./p0.^3 ...
            - (E0E+p.^2).*Eps./p.^3 + 2*k.*E0E./p0p.^2 );

        d1sigma(i) = PreFactor .* ( Term1 + Term2 + Term3 ); %.*(k<(E0-1));
    end
        spec =  (dens*(Z.^2)')*d1sigma;
end

%function pow = power(nspecies, Z, density, r0) %Needed for testung 4BS,
%however for some unknown reason, does not work, and causes other functions
%to also not work. It gets called when it shouldn't
%    Ip = zeros(1, length(nspecies));
%    k = 1; ip = 1;
%    for jp=1:length(nspecies)
%        while ip <= nspecies(jp) + k-1
%            Z(ip)
%            Ip(jp) = Ip(jp) + density(ip)*Get_4BS(Z(ip), r0);
%            ip = ip +1;
%        end
%        k = ip;
%    end
%    pow = Ip;
%end
close all;
figure;
semilogy(r(end,:),'LineWidth',2);hold on;
semilogy(R(:),'LineWidth',2);

figure;

for j=1:N_gen
    
    scatter(t,r(:,j)/(R0^3/(2^(j-1)))^(1/3),'.'); hold on;
    
end

xlabel time(days)
ylabel R/R0

phie = Me'./(Me(:)+Mc(:)+Mm(:));
phim = Mm'./(Me(:)+Mc(:)+Mm(:));
phic = Mc'./(Me(:)+Mc(:)+Mm(:));

figure;
plot(phie,'LineWidth',2);hold on;
plot(phim,'LineWidth',2);
plot(phic,'LineWidth',2);
xlabel generation
ylabel \phi
legend('elastin','SMC','Collagen');

figure;
hr = h./r(end,:);
stress = p_mid(:)./hr';

plot(hr,'LineWidth',2);
xlabel generation
ylabel h/r
figure;
plot(stress,'LineWidth',2);
ylabel stress(Pa)
xlabel generation

figure;
plot(tau,'LineWidth',2)
ylabel WSS(Pa)
xlabel generation

for k=1:N_gen-1
ksi(k) = 1/(log2(r(end-2,k))-log2(r(end-2,k+1)));
end
figure;
plot(ksi,'LineWidth',2);
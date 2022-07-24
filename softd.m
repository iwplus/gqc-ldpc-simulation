################################################################################
##### program utama untuk Log-Likelihood Ratios Sum-Product Algorithm ##########
################################################################################

### (1) Inisialisasi matriks H ####

partisi = randi([2 4],1,3); #permutasi dari kode dimana partisi paling ujung sama dengan ko-dimensi
disp('Cycle structure = ');
disp(partisi);
disp('');
disp('');
n = sum(partisi); #panjang kode
disp('Panjang kode = ');
disp(n);
disp('');
disp('');
g = randi([0 1], 1, n); #baris pertama matriks H
r = partisi(length(partisi)); #ko-dimensi kode LDPC
disp('Ko-dimensi kode = ');
disp(r);
disp('');
disp('');
H=generateh(partisi,g,r); #matriks H
disp('H = ');
disp(H);
disp('');
disp('');

###### (2) Enkoding ####
#### catatan : pastikan blok paling kanan matriks H sudah tak-singular ####

### Hitung matriks G ####

indawal = sum(partisi)-partisi(length(partisi))+1;
indakhir = indawal+partisi(length(partisi))-1;
Hk = H(:,indawal:indakhir);
Hkinv = g2rref(Hk);
G = eye(n-r);
ind = 0;
Gtemp=[];
for i = 1:(length(partisi)-1),
  if i == 1,
    B1 = Hkinv*H(1:r,1:partisi(i));
  else,
    sind = sum(partisi(1:(i-1)));
    ind1 = sind+1;
    ind2 = ind1+partisi(i)-1;
    B1 = Hkinv*H(1:r,ind1:ind2);
  endif
  B2 = mod(B1,2);
  B3 = transpose(B2);  
  Gtemp = [Gtemp;B3];
endfor

G = [G Gtemp];
disp('G = ');
disp(G);
disp('');
disp('');

#### generate pesan asli secara random ####
m = zeros(1,n-r);
nsatu = randperm(n-r,1);
indsatu = randperm(n-r,nsatu);

for iterisim = 1:length(indsatu),
  m(indsatu(iterisim)) = 1;
endfor
##disp('Pesan asli = ');
##disp(m);
##disp('');
##disp('');

x = mod(m*G,2); #### hasil enkoding
##disp('pesan yang dikirim = ');
##disp(x);
##disp('');
##disp('');

#t = 1; #### bobot vektor error
ber = [];
numerror = [];
bataserror = ceil(2*n/3);
for t = 1:bataserror, 
inderror = randperm(n,t);
xerror = x;

for i = 1:t, #### penambahan error ke pesan yang dikirim
  xerror(inderror(i)) = mod(xerror(inderror(i))+1,2);
endfor

##disp('pesan yang diterima = ');
##disp(xerror);
##disp('');
##disp('');

### (3) Inisialisasi LLR-SPA ###

y = xerror;
Gamma = zeros(n,r);
Lambda = zeros(r,n);

for i = 1:n,
  for j = 1:r,
    if H(j,i) == 1,
      Gamma(i,j) = llr(y(i),t,n);
    endif
  endfor
endfor

##disp('Gamma awal = ');
##disp(Gamma);
##disp('');
##disp('');

iterasimaks = 2;
iter = 0;
while iter <= iterasimaks, #### iterasi LLR-SPA
  
#### (4) Left semi-iteration #####

for k = 1:r,
  Ak = find(H(k,:)); #### A(k) = daftar semua tetangga dari check node c_k
  
  for i = 1:length(Ak),
    Aki = Ak;
    Aki(i) = [];
    tanhtemp = 1;
    for j = 1:length(Aki),
      tanhtemp = tanhtemp*tanh(1/2*Gamma(Aki(j),k));
    endfor
    Lambda(k,Ak(i)) = 2*atanh(tanhtemp);
  endfor
endfor

##disp('Lambda terupdate = ');
##disp(Lambda);
##disp('');
##disp('');

#### (5) Right semi-iteration #####
Gammadec = zeros(1,n);
for i = 1:n,
  Bi = transpose(find(H(:,i))); ###B(i) = daftar semua tetangga dari var node v_i
  
  for k = 1:length(Bi),
    Bik = Bi;
    Bik(k) = [];
    Gamtemp = 0;
    for j = 1:length(Bik),
      Gamtemp = Gamtemp + Lambda(Bik(j),i);
    endfor
    Gamma(i,Bi(k)) = llr(y(i),t,n) + Gamtemp;
  endfor
  Gamtemp2 = 0;
  for j = 1:length(Bi),
      Gamtemp2 = Gamtemp2 + Lambda(Bi(j),i);
  endfor 
  Gammadec(i) = llr(y(i),t,n) + Gamtemp2;  
endfor

##disp('Gamma terupdate = ');
##disp(Gamma);
##disp('');
##disp('');

#disp('Gamma_i untuk keputusan = ');
#disp(Gammadec);
#disp('');
#disp('');


######### (4) Decision ##########

xcorrect = zeros(1,n);
for i = 1:n,
  if Gammadec(i) >= 0,
    xcorrect(i) = 0;
  else,
    xcorrect(i) = 1;
  endif
endfor

#disp('x hasil koreksi = ');
#disp(xcorrect);
#disp('');
#disp('');

if xcorrect == x,
  iter = iterasimaks + 10;
else,
  iter = iter + 1;
endif  

endwhile

##disp('x hasil koreksi = ');
##disp(xcorrect);
##disp('');
##disp('');

######### hitung Bit Error Rate (BER) ####

selisihv = mod(x-xcorrect,2);
selisihind = find(selisihv);
ber = [ber length(selisihind)/n];
numerror = [numerror t];

endfor

###### plot BER terhadap banyaknya error yang ditambahkan ######

plot(numerror,ber,'-s', 'MarkerSize',10,'MarkerFaceColor',[1 .6 .6])
set(gca,'Ylim',[0 1])
set(gca,'Xlim',[0 n])
xlabel('Number of errors')
ylabel('Bit Error Rate')
title('Bit Error Rate vs Error added')
grid on
%this program computes approximate output distributions from boson
%samplers, and generates samples from those distributions using a monte
%carlo method. If you use this code, please cite Renema et al. 
%Phys. Rev. Lett. 120 (22), 220502 and arXiv:1809.01953, which are the
%underlying works.
clear;
clc;

%specify unitary matrix
U = randU(10);

%specify input modes
inputmodes = [1 2 3 4];

%specify mutual distinguishability & transmission
eta = 0.5;
x = 0.8;

%specify truncation
kmax = 2;

%specify number of output photons (not sampling over those for the moment)
k = 3;

%control parameters: sample from output distribution, compute output
%distribution. Give number of samples
wantdist = 1;
wantsample = 1;
Nsample = 1e5;

%Here begins the actual computation. Do not modify beyond this point unless 
%you know what you are doing.
n = size(inputmodes,2);
N = size(U,1)

Nlist = 1:N;
Nlist = Nlist';

%generate list of all no-collision outputs
for ii = 1:k-1
    newlist = [];
    for jj = 1:size(Nlist,1)
        for kk = 1:N %Nlist(jj,end)
            newlist(end+1,:) = [Nlist(jj,:) kk];
        end
    end
    Nlist = newlist
end

Pexact = zeros(size(Nlist,1),1);
Papprox = zeros(size(Nlist,1),1);

%generate all possible loss mechanisms
combs = nchoosek(inputmodes,k)

if wantdist
    for ii = 1:size(Nlist,1)
        for jj = 1:size(combs,1)
            perm = perms(combs(jj,:));
            for kk = 1:size(perm,1)
                korder = k - sum(perm(kk,:) == combs(jj,:));
                pterm = x^korder * permanentRyser(U(combs(jj,:),Nlist(ii,:)).*conj(U(perm(kk,:),Nlist(ii,:))))/(nchoosek(n,k)*factorial(k));
                Pexact(ii) = Pexact(ii) + pterm;
                if korder < kmax+1
                    Papprox(ii) = Papprox(ii) + pterm;
                end
            end
        end
    end

    %remove imaginary part of Papprox, Pexact, which arises due to rounding errors
    Papprox = real(Papprox);
    Pexact = real(Pexact);
    
    figure(1)
    plot([Papprox Pexact])
 
end

%now sample over the distribution
if wantsample
    b = Nlist(1,:)
    comb = inputmodes(1:k)
    proba = Pexact(1);
    currentpos = 1;
    lsample = [];

    %slist = [];
    bprob = zeros(size(Nlist,1),1);

    for ii = 1:Nsample
        if ii/1000 == floor(ii/1000)
            ii
        end

        testpos = randi(size(Nlist,1));
        btest = Nlist(testpos,:);
        combtest = combs(randi(size(combs,1)),:);
        perm = perms(combtest);
        testprob = 0;

        for jj = 1:size(perm,1)
            if (k - sum(perm(jj,:) == combtest)) < kmax+1
                testprob = testprob + permanentRyser(U(combtest,btest).*conj(U(perm(jj,:),btest)));
            end
        end

        if testprob > proba
            b = btest;
            proba = testprob;
            comb = combtest;
            currentpos = testpos;
        else
            if testprob / proba > rand(1)
                b = btest;
                proba = testprob;
                comb = combtest;
                currentpos = testpos;
            end
        end
        bprob(currentpos) = bprob(currentpos) + 1/Nsample;
        lsample(end+1) = currentpos;
    end

        figure(2)
        plot([Papprox Pexact bprob])
end
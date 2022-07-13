function [k,lam,S_mu,S_a,EmptyElList] = param_decoup_main(F_t,og_matprop)
    h = waitbar(0,'Progress: 0%');
    for i = 1:length(F_t)
        [k{i},lam{i},S_mu{i},S_a{i},EmptyElList{i}] = param_decoup(F_t{i},og_matprop);
        waitbar(i/length(F_t),h,...
            ['Progress for run #',num2str(i),'/',num2str(length(F_t)),': ',num2str(floor(100*i/length(F_t))),'%'])
    end
    close(h)
end

function [k,lam,S_mu,S_a,EmptyElList] = param_decoup(F_t,og_matprop)
    digits(4);
    sz = size(F_t{1,1});
    if sz(1) == 1
        def_grad = convJonFtoDennisF_lin(F_t);
    else
        def_grad = convJonFtoDennisF(F_t);
    end
    def_grad2 = def_grad(:);
    ppm = ParforProgressbar(numel(def_grad));
    EmptyElList = [];
    mu = og_matprop(1,:);
    a = og_matprop(2,:);
    lam = NaN(size(def_grad2));
    k = NaN(size(def_grad2));
    S_mu_1 = NaN(size(def_grad2));
    S_a_1 = NaN(size(def_grad2));
    S_mu_2 = NaN(size(def_grad2));
    S_a_2 = NaN(size(def_grad2));
    S_mu_3 = NaN(size(def_grad2));
    S_a_3 = NaN(size(def_grad2));
    S_mu_4 = NaN(size(def_grad2));
    S_a_4 = NaN(size(def_grad2));
    parfor ii = 1:numel(def_grad)
        if sum(isnan(def_grad2{ii})) == 0
            F = def_grad2{ii};
            F = F/(det(F)^(1/3));
            [lam_temp,k_temp] = FToLamAndK_v2(F);
            if length(lam_temp) == 0 || isreal(lam_temp) == 0 || isreal(k_temp) == 0
                EmptyElList = [EmptyElList ii];
            else
                lam(ii) = lam_temp; k(ii) = k_temp;
                [S_mu_1(ii),S_a_1(ii)] = sensitivity(k(ii),lam(ii),a(1),mu(1));
                [S_mu_2(ii),S_a_2(ii)] = sensitivity(k(ii),lam(ii),a(2),mu(2));
                [S_mu_3(ii),S_a_3(ii)] = sensitivity(k(ii),lam(ii),a(3),mu(3));
                [S_mu_4(ii),S_a_4(ii)] = sensitivity(k(ii),lam(ii),a(4),mu(4));
            end
        end
        ppm.increment();
    end
    if sz(1) ~= 1
        % Experimental
        lam = convDennisStoJonS(reshape(lam,size(def_grad)));
        k = convDennisStoJonS(reshape(k,size(def_grad)));
        S_mu = {S_mu_1,S_mu_2,S_mu_3,S_mu_4};
        S_a = {S_a_1,S_a_2,S_a_3,S_a_4};
        for i = 1:length(a)
            S_mu{i} = convDennisStoJonS(reshape(S_mu{i},size(def_grad)));
            S_a{i} = convDennisStoJonS(reshape(S_a{i},size(def_grad)));
        end
    else
        % Simulations
        lam = convDennisStoJonS_lin(reshape(lam,size(def_grad)));
        k = convDennisStoJonS_lin(reshape(k,size(def_grad)));
        S_mu = {S_mu_1,S_mu_2,S_mu_3,S_mu_4};
        S_a = {S_a_1,S_a_2,S_a_3,S_a_4};
        for i = 1:length(a)
            S_mu{i} = convDennisStoJonS_lin(reshape(S_mu{i},size(def_grad)));
            S_a{i} = convDennisStoJonS_lin(reshape(S_a{i},size(def_grad)));
        end
    end
    delete(ppm)
end

function [S_mu, S_a] = sensitivity(k,lam,a,mu)
    S_mu = (-k^2+1.5*k+0.5)*lam^(a*(-k^2 ...
        +1.5*k+0.5)-1) + (-k+0.5)*lam^(a*(-k+0.5)-1) + ...
        (k^2-0.5*k-1)*lam^(a*(k^2-0.5*k-1)-1);
    S_a = mu*(log(lam)/lam)*...
        ((-k^2+1.5*k+0.5)^2)*lam^(a*(-k^2 ...
        +1.5*k+0.5)) + ((-k+0.5)^2)*lam^(a*(-k+0.5)) + ...
        ((k^2-0.5*k-1)^2)*lam^(a*(k^2-0.5*k-1));
end
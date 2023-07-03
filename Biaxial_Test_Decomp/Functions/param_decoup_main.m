function [k,lam,EmptyElList] = param_decoup_main(F_t,parallel,decomp_mode)
    for i = 1:length(F_t)
        [k{i},lam{i},EmptyElList{i}] = param_decoup(F_t{i},parallel,decomp_mode);
    end
end

function [k,lam,EmptyElList] = param_decoup(F_t,parallel,decomp_mode)
    digits(4);
    sz = size(F_t{1,1});
    if sz(1) == 1
        def_grad = convJonFtoDennisF_lin(F_t);
    else
        def_grad = convJonFtoDennisF(F_t);
    end
    def_grad2 = def_grad(:);
    EmptyElList = [];
    lam = NaN(size(def_grad2));
    k = NaN(size(def_grad2));
    switch parallel
        case 1
            parfor ii = 1:numel(def_grad)
                if sum(isnan(def_grad2{ii})) == 0
                    F = def_grad2{ii};
                    F = F/(det(F)^(1/3));
                    switch decomp_mode
                        case 'ftolamandk'
                            [lam_temp,k_temp] = FToLamAndK_v2(F);
                        case 'ftok2andk3'
                            [~,lam_temp,k_temp] = FToK2AndK3(F);
                    end
                    if isempty(lam_temp) == true || isreal(lam_temp) == 0 || isreal(k_temp) == 0
                        EmptyElList = [EmptyElList ii];
                    else
                        lam(ii) = lam_temp; k(ii) = k_temp;
                    end
                end
            end
        case 0
            for ii = 1:numel(def_grad)
                if sum(isnan(def_grad2{ii})) == 0
                    F = def_grad2{ii};
                    F = F/(det(F)^(1/3));
                    switch decomp_mode
                        case 'ftolamandk'
                            [lam_temp,k_temp] = FToLamAndK_v2(F);
                        case 'ftok2andk3'
                            [~,lam_temp,k_temp] = FToK2AndK3(F);
                    end
                    if isempty(lam_temp) == true || isreal(lam_temp) == 0 || isreal(k_temp) == 0
                        EmptyElList = [EmptyElList ii];
                    else
                        lam(ii) = lam_temp; k(ii) = k_temp;
                    end
                end
            end
    end
    if sz(1) ~= 1
        % Experimental
        lam = convDennisStoJonS(reshape(lam,size(def_grad)));
        k = convDennisStoJonS(reshape(k,size(def_grad)));
    else
        % Simulations
        lam = convDennisStoJonS_lin(reshape(lam,size(def_grad)));
        k = convDennisStoJonS_lin(reshape(k,size(def_grad)));
    end
end
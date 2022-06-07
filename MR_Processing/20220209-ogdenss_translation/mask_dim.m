load("mask_data_850.mat") % mask, mask_, maskthresh, and osc from runComplex
units = 'in';
% units = 'mm';

[s1,s2,s3] = size(mask_);

ct_1 = 1; ct_2 = 1;
for i = 1:s2
    for j = 1:s3
        thick_samp_1 = sum(mask_(1:s1/2,i,j));
        thick_samp_2 = sum(mask_(s1/2:s1,i,j));
        if thick_samp_1 > 0
            samp_1{1}(ct_1) = thick_samp_1;
            ct_1 = ct_1 + 1;
        end
        if thick_samp_2 > 0
            samp_2{1}(ct_2) = thick_samp_2;
            ct_2 = ct_2 + 1;
        end
    end
end
figure
histogram(samp_1{1})
xlim([0,max([samp_1{1} samp_2{1}])])
title('Sample Thickness')
xlabel('Voxel length')
hold on
histogram(samp_2{1})
legend('Sample 1','Sample 2')

cw_1 = 1; cw_2 = 1;
for i = 1:s1
    for j = 1:s3
        if i <= s1/2
            width_samp_1 = sum(mask_(i,:,j));
            if width_samp_1 > 0
                samp_1{2}(cw_1) = width_samp_1;
                cw_1 = cw_1 + 1;
            end
        else
            width_samp_2 = sum(mask_(i,:,j));
            if width_samp_2 > 0
                samp_2{2}(cw_2) = width_samp_2;
                cw_2 = cw_2 + 1;
            end
        end
    end
end
figure
histogram(samp_1{2})
xlim([0,max([samp_1{2} samp_2{2}])])
title('Sample Width')
xlabel('Voxel length')
hold on
histogram(samp_2{2})
legend('Sample 1','Sample 2')

cl_1 = 1; cl_2 = 1;
for i = 1:s1
    for j = 1:s2
        if i <= s1/2
            length_samp_1 = sum(mask_(i,j,:));
            if length_samp_1 > 0
                samp_1{3}(cl_1) = length_samp_1;
                cl_1 = cl_1 + 1;
            end
        else
            length_samp_2 = sum(mask_(i,j,:));
            if length_samp_2 > 0
                samp_2{3}(cl_2) = length_samp_2;
                cl_2 = cl_2 + 1;
            end
        end
    end
end
figure
histogram(samp_1{3})
xlim([0,max([samp_1{3} samp_2{3}])])
title('Sample Length')
xlabel('Voxel length')
hold on
histogram(samp_2{3})
legend('Sample 1','Sample 2')

samp_1{4} = [mean(samp_1{1})*osc(1) mean(samp_1{2})*osc(2) mean(samp_1{3})*osc(3)];
samp_2{4} = [mean(samp_2{1})*osc(1) mean(samp_2{2})*osc(2) mean(samp_2{3})*osc(3)];

samp_1{5} = [max(samp_1{1})*osc(1) max(samp_1{2})*osc(2) max(samp_1{3})*osc(3)];
samp_2{5} = [max(samp_2{1})*osc(1) max(samp_2{2})*osc(2) max(samp_2{3})*osc(3)];

samp_1{6} = [std(samp_1{1})*osc(1) std(samp_1{2})*osc(2) std(samp_1{3})*osc(3)];
samp_2{6} = [std(samp_2{1})*osc(1) std(samp_2{2})*osc(2) std(samp_2{3})*osc(3)];

if units == 'in'
    samp_1{4} = samp_1{4}/25.4;
    samp_2{4} = samp_2{4}/25.4;
    samp_1{5} = samp_1{5}/25.4;
    samp_2{5} = samp_2{5}/25.4;
    samp_1{6} = samp_1{6}/25.4;
    samp_2{6} = samp_2{6}/25.4;
end

Property = {'Mean';'Max';'Std'};
A1 = [samp_1{4};samp_1{5};samp_1{6}];
Thickness = A1(:,1); Width = A1(:,2); Length = A1(:,3);

Sample1 = table(Property,Thickness,Width,Length)

A2 = [samp_2{4};samp_2{5};samp_2{6}];
Thickness = A2(:,1); Width = A2(:,2); Length = A2(:,3);

Sample2 = table(Property,Thickness,Width,Length)
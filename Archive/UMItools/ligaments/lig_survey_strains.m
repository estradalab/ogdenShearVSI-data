function h=lig_survey_strains(varargin)
filename=varargin{1};
%%
P=load(filename);  %use this one (5mm)
L=P.Lagrange;

%%
hf=0;

mask=P.mask_hr.*P.mask; %and the masks of the GE3D data and the DESTE data

%identify a default center of the sample, using the mask;
si=size(mask);
mc1=round(si(1)/2);
testslice=squeeze(mask(mc1,:,:));
mc2=round(sum(testslice,2)'*(1:si(2))'/sum(testslice(:)));
mc3=round((1:si(3))*sum(testslice,1)'/sum(testslice(:)));

mask(mask==0)=NaN;
cmap=jet(256);


%% top view Lagrangian strains
for sli=mc3+(-21:3:21);
    Lagrangian_derivatives_flag=true;
    if Lagrangian_derivatives_flag;
        hf=hf+1;
        h(hf)=figure('position',[10 10 1900 950]);
        set(gcf,'PaperPositionMode','manual','colormap',cmap,'PaperOrientation','landscape');
        for enc=1:3;
            tsubplot(3,5,1+(enc-1)*5,8);
            imagesc(P.axis2,P.axis1,angle(P.data(:,:,sli,enc)).*mask(:,:,sli),[-pi pi]);
            axis image
            grid on
            set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-');%,'Yticklabel','','Xticklabel','');
            %ylabel(['slice ' num2str(sli)]);
            colorbar;
            switch enc
                case 1
                    title(['phase (ud)']);
                case 2
                    title(['phase (lr)']);
                case 3
                    title(['phase (io)']);
            end
            
            
            for jj=1:3;
                tsubplot(3,5,jj+1+(enc-1)*5,8);
                %imagesc(P.axis2,P.axis1,blur(P.strains(:,:,sli,enc,jj)),[-0.2 0.2]);
                imagesc(P.axis2,P.axis1,(L(:,:,sli,enc,jj)).*mask(:,:,sli),[-0.20 0.20]);
                axis image
                grid on
                set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
                title(['L ' num2str(enc) num2str(jj)],'color','k');
                hc=colorbar;
%                 pv=get(hc,'Position');
%                 pv(3)=0.005;
%                 set(hc,'Position',pv);
            end
            
            
            %plot an axial cut, with a line indicating the slice
            tsubplot(3,5,15,8);
            % plot the magnitudes in the 5th column
            MM=permute(squeeze(log10(P.HIRES.magnitude(mc1,:,:)).*mask(mc1,:,:)),[2 1]);
            imagesc(P.HIRES.axis2,P.HIRES.axis3,MM);
            axis image
            set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','YDir','normal');%,'Yticklabel','','Xticklabel',''
            title(['log(|HIRES|), slide#' num2str(hf)]);
            hold on;
            plot([P.HIRES.axis2(1) P.HIRES.axis2(end)],P.HIRES.axis3(sli)*[1 1],'w'); 
            xlabel('dir2 (mm)');
            ylabel('dir3 (mm)');
            colorbar;
            
            
            % plot the magnitudes in the 5th column
            tsubplot(3,5,5,8);
            imagesc(P.HIRES.axis2,P.HIRES.axis1,log10(P.HIRES.magnitude(:,:,sli)).*mask(:,:,sli));
            axis image
            grid on
            set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
            title('log(|HIRES|)');
            colorbar;
            
            % plot the magnitudes in the 5th column
            tsubplot(3,5,10,8);
            imagesc(P.axis2,P.axis1,log10(abs(P.data(:,:,sli,1))).*mask(:,:,sli));
            axis image
            grid on
            set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','');
            title('log(|DESTE|)');
            colorbar;
            orient landscape            
            
        end
    end
end

%% cross view Lagrangian strains
for sli=mc1+(-75:5:75);
    Lagrangian_derivatives_flag=true;
    if Lagrangian_derivatives_flag;
        hf=hf+1;
        h(hf)=figure('position',[50 10 1900 950]);
        set(gcf,'PaperPositionMode','manual','colormap',cmap,'PaperOrientation','landscape');
        for enc=1:3;
            tsubplot(3,5,1+(enc-1)*5,8);
            imagesc(P.axis2,P.axis3,permute(squeeze(angle(P.data(sli,:,:,enc)).*mask(sli,:,:)),[2 1]),[-pi pi]);
            axis image
            grid on
            set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','YDir','normal');%,'Yticklabel','','Xticklabel',''
            ylabel(['slice ' num2str(sli)]);
            hc=colorbar('Location','SouthOutside');
%             pv=get(hc,'Position');
%             pv(3)=0.005;
%             set(hc,'Position',pv);
%             
            switch enc
                case 1
                    title(['phase 1 (io)']);
                case 2
                    title(['phase 2 (lr)']);
                case 3
                    title(['phase 3 (ud)']);
            end
            
            
            for jj=1:3;
                tsubplot(3,5,jj+1+(enc-1)*5,8);
                imagesc(P.axis2,P.axis3,permute(squeeze((L(sli,:,:,enc,jj)).*mask(sli,:,:)),[2 1]),[-0.2 0.2]);
                axis image
                grid on
                set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','','YDir','normal');
                title(['L ' num2str(enc) num2str(jj)]);
                hc=colorbar('Location','SouthOutside');
%                 pv=get(hc,'Position');
%                 pv(3)=0.005;
%                 set(hc,'Position',pv);
            end
            
            % plot the magnitudes in the 5th column
            MM=permute(squeeze(log10(P.HIRES.magnitude(sli,:,:)).*mask(sli,:,:)),[2 1]);
            tsubplot(3,5,5,8);
            imagesc(P.HIRES.axis2,P.HIRES.axis3,MM);
            axis image
            grid on
            set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','','YDir','normal');
            title('log(|HIRES|)');
            colorbar('Location','SouthOutside');
            
            % plot the magnitudes in the 5th column
            MM=permute(squeeze(log10(abs(P.data(sli,:,:,1))).*mask(sli,:,:)),[2 1]);
            tsubplot(3,5,10,8);
            imagesc(P.axis2,P.axis3,MM);
            axis image
            grid on
            set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-','Yticklabel','','Xticklabel','','YDir','normal');
            title('log(|DESTE|)');
            colorbar('Location','SouthOutside');
            
            % plot the magnitudes in the 5th column
            tsubplot(3,5,15,8);
            imagesc(P.HIRES.axis2,P.HIRES.axis1,log10(P.HIRES.magnitude(:,:,mc3)).*mask(:,:,mc3));
            axis image
            set(gca, 'XColor', 'k','Ycolor','k','gridlinestyle','-');%,'Yticklabel','','Xticklabel','');
            title(['log(|HIRES|), slide#' num2str(hf)]);
            hold on;
            plot([P.HIRES.axis2(1) P.HIRES.axis2(end)], P.HIRES.axis1(sli)*[1 1],'w');
            xlabel('dir2 (mm)');
            ylabel('dir1 (mm)');
            colorbar;
            
            
        end
    end
end
if any(strcmp(varargin,'pdf'));
    pdfappend(h,[filename '.pdf'],'size',[10 8]);
end

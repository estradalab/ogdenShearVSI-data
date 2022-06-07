%test out resolution vs lambda, all in the RO direction

%lambda=1;

clear;


% timestampvec={'1254','1303','1312','1321','1330','1339','1350','1401','1410',...
%     '1419','1428','1435','1448','1456','1504','1512','1520','1526','1533',...
%     '1539'};


%files for lamda over rox variability
%timestampvec={'1134','1142','1150','1158','1206','1214','1222','1230','1237'};

%files for resolution variability, upsampled to common 256 readout
timestampvec={'1448','1456','1350','1504','1512','1520','1526','1006', '1533','1539'};

for jj=1:numel(timestampvec);
    dummy=lig_deste_1d(timestampvec{jj},'pixvec',[256 32 20],'apod');
%     data{jj}=dummy.origdata;
%     
%     linear approximation of the phase wind, using a central patch of  the
%     central slice:
%     test=data{jj}(50:90,7:13,13)./data{jj}(51:91,7:13,13);
%     slope=mean(angle(test(:)));
%     
%     si=size(data{jj});
%     correction=exp(1i*(1:si(1))'*slope * ones(1,si(2)));
%     
%     test2=data{jj}(50:90,7:13,13).*correction(50:90,7:13);
%     offset=mean(angle(test2(:)));
%     
%     correction=correction*exp(-1i*offset);
%     
%     
%     
%     figure('Position',[100 100 2400 1200],'name',dummy.comment);
%     for jjj=1:8;
%         
%         slice_index=jjj+7;
%         
%         subplot(3,10,1+(jjj-1));
%         imagesc(angle(data{jj}(:,:,slice_index)));
%         title(['slice ' num2str(slice_index)]);
%         
%         subplot(3,10,11+(jjj-1));
%         imagesc(angle(data{jj}(:,:,slice_index).*correction));
%         title(['slice ' num2str(slice_index)]);
%         
%         
%         
%         subplot(3,10,21+(jjj-1));
%         %imagesc(log(abs(data{jj})));
%         imagesc((abs(data{jj}(:,:,slice_index))),[0 7000]);
%         title(['slice ' num2str(slice_index)])
%     end
%     
%     meanstrain=slope*dummy.lambda/2/pi/abs(dummy.axis1(1)-dummy.axis1(2));
%     
%     subplot(3,5,5);
%     plot(squeeze(angle(data{jj}(:,8,11:14))));
%     title('phase');
%     
%     subplot(3,5,10);
%     plot(angle(squeeze(data{jj}(:,8,11:14)).*(correction(:,8)*ones(1,4))));
%     title(['phase minus ' num2str(round(slope/pi*180)) ' deg/pix, strain=' num2str(round(meanstrain*1000)/1000)]);
%     
%     subplot(3,5,15);
%     plot(squeeze(abs(data{jj}(:,8,11:14))));
%     set(gca,'ylim',[0 7000]);
%     title('magnitude')
%     pause(1);
end
function varargout = colormanifolds(play)
% Different ways of visualising how the combining of colours (i.e.
% modifying a colour in two dimensions) translates to its effects in 3
% dimensional colour space (either HSV, xyY, etc.). This function requires
% the toolbox optprop by Jerker Wagberg.
%
% Joshua Harvey 2017

if nargin == 0
    play = 'general';
end

v = 0;

if strcmp(play,'general')
    %% General implementation in HSV space
    
    HH1=[0 1/3 2/3 1/6 3/6 5/6];
    HH2=[1/3 2/3 0 3/6 5/6 1/6];
    sats=fliplr([1 .9 .75 .5 .25 .1 .05]);
    
    % % just for red-green
    % HH1=[0];
    % HH2=[1/3];
    % sats=[1];
    
    % % wide range of all colours
    % HH1 = [.1 .3 .5 .7 .9 .1 .3 .5 .6 .8];
    % HH2 = [HH1(3:7) HH1(3:7)-.1];
    % sats = [1 .9 .8 .7 .6 .3 .15 0];
    
    % % monohue (fully saturated)
    % HH1 = linspace(0,.99,100);
    % HH2 = HH1+.005;
    % sats = 1;
    
    % % oponent colours
    % HH1 = linspace(0,1,10);
    % HH2plus = HH1+.5;
    % integ=floor(HH2plus);
    % HH2=HH2plus-integ;
    % sats = 1;
    
    sscale = .5*pi; % .5*pi for cylinder
    
    figure('Color',[0 0 0]);
    v = 0;
    for j = 1:length(sats)
        for i = 1:length(HH1)
            % define colors in HSV
            H1 = hsv2rgb([HH1(i),sats(j),1]); % red
            H2 = hsv2rgb([HH2(i),sats(j),1]); % yellow
            
            h1r = linspace(0,H1(1),100);
            h1g = linspace(0,H1(2),100);
            h1b = linspace(0,H1(3),100);
            h2r = linspace(0,H2(1),100);
            h2g = linspace(0,H2(2),100);
            h2b = linspace(0,H2(3),100);
            
            [RX,RY] = meshgrid(h1r,h2r);
            [GX,GY] = meshgrid(h1g,h2g);
            [BX,BY] = meshgrid(h1b,h2b);
            
            R = (RX + RY); % add components of colors
            G = (GX + GY);
            B = (BX + BY);
            
            RGB = cat(3,R,G,B);
            RGBhsv = rgb2hsv(RGB);
            [HSVX, HSVY, HSVZ] = sph2cart(RGBhsv(:,:,1)*2*pi,(pi/2)-RGBhsv(:,:,2),0.5*RGBhsv(:,:,3));
            
            % need to correct V
            
            RGB1 = cat(3,RX,GX,BX); RGB2 = cat(3,RY,GY,BY);
            HSV1 = rgb2hsv(RGB1); HSV2 = rgb2hsv(RGB2);
            V12 = .5*(HSV1(:,:,3)+HSV2(:,:,3));
            [HSVX, HSVY, HSVZ] = sph2cart(RGBhsv(:,:,1)*2*pi,(pi/2)-sscale*RGBhsv(:,:,2),V12);
            
            % color displayed
            mapper = rgb2hsv(RGB/2);
            mapper2 = mapper;
            mapper2(:,:,3) = mapper(:,:,3)*2;
            RGBmapper = hsv2rgb(mapper2);
            
            mapper = hsv2rgb(cat(3,RGBhsv(:,:,1:2),V12));
            
            %         % rescale radius by value
            %         HSVZ = RGBhsv(:,:,3);
            
            s=surf(HSVX,HSVY,HSVZ,mapper); shading flat
            s.FaceAlpha=j/10+.3;
            xlim([-1 1]);
            ylim([-1 1]);
            zlim([-.2 1]);
            hold on
            axis off
            view([v 30])
            pause(0.1)
            v=v+5;
            ax=gca;
            ax.DataAspectRatio = [2 2 1];
            axis vis3d;
            if i*j == 1
                camzoom(1.3)
            end
        end
    end % for HSV space
    
elseif strcmp(play,'surface')
    %% define colors in hsv and plot in two dimensional surface
    r = [0 1 1];
    g = [1/3 1 1];
    o = ones([1,100]);
    
    rs = repmat(r,[100 1]).*[o;o;linspace(0,1,100)]';
    gs = repmat(g,[100 1]).*[o;o;linspace(0,1,100)]';
    
    [Hr,Hg] = meshgrid(rs(:,1),gs(:,1));
    [Sr,Sg] = meshgrid(rs(:,2),gs(:,2));
    [Vr,Vg] = meshgrid(rs(:,3),gs(:,3));
    
    rgbRG = hsv2rgb(Hr,Sr,Vr) + hsv2rgb(Hg,Sg,Vg);
    hsvRG = rgb2hsv(rgbRG);
    hueRG = hsvRG(:,:,1);
    
    surf(hueRG,rgbRG)
    shading flat
    
    hold on
    
    % add another colour
    g = [1/3 1 1];
    r = [2/3 1 1];
    o = ones([1,100]);
    
    rs = repmat(r,[100 1]).*[o;o;linspace(0,1,100)]';
    gs = repmat(g,[100 1]).*[o;o;linspace(0,1,100)]';
    
    [Hr,Hg] = meshgrid(rs(:,1),gs(:,1));
    [Sr,Sg] = meshgrid(rs(:,2),gs(:,2));
    [Vr,Vg] = meshgrid(rs(:,3),gs(:,3));
    
    rgbRG = hsv2rgb(Hr,Sr,Vr) + hsv2rgb(Hg,Sg,Vg);
    hsvRG = rgb2hsv(rgbRG);
    hueRG = hsvRG(:,:,1);
    figure('Color',[0 0 0]);
    surf(hueRG,rgbRG)
    shading flat
    
elseif strcmp(play,'physical')
    %% physical colours (wavelengths) %490-700 are opposite
    
    % %nearly monochromatic
    %waves1 = linspace(380,700,100);
    %waves2 = waves1+3;
    waves = linspace(380,700,40);
    waves1 = [waves waves waves waves waves];
    waves2 = rem([waves+5,waves+20,waves+40,waves+150,waves+210],700);
    
    sscale = 1; %.5*pi; % .5*pi for cylinder
    
    figure('Color',[0 0 0]);
    a = 148;
    e = 0;
    for j = 1 %:length(sats)
        for i = 1:length(waves1)
            % define colors in RGB
            H1 = MakeBruntonsRGB(waves1(i)); %
            H2 = MakeBruntonsRGB(waves2(i)); %
            
            h1r = linspace(0,H1(1),100);
            h1g = linspace(0,H1(2),100);
            h1b = linspace(0,H1(3),100);
            h2r = linspace(0,H2(1),100);
            h2g = linspace(0,H2(2),100);
            h2b = linspace(0,H2(3),100);
            
            
            [RX,RY] = meshgrid(h1r,h2r);
            [GX,GY] = meshgrid(h1g,h2g);
            [BX,BY] = meshgrid(h1b,h2b);
            
            R = (RX + RY); % add components of colors
            G = (GX + GY);
            B = (BX + BY);
            
            RGB = cat(3,R,G,B); % imshow(RGB)
            RGBhsv = rgb2hsv(RGB);
            % RGBhsv(:,:,2) = ones(100); % set all sats to 1
            [HSVX, HSVY, HSVZ] = sph2cart(RGBhsv(:,:,1)*2*pi,(pi/2)-RGBhsv(:,:,2),0.5*RGBhsv(:,:,3));
            
            % need to correct V
            
            RGB1 = cat(3,RX,GX,BX); RGB2 = cat(3,RY,GY,BY);
            HSV1 = rgb2hsv(RGB1); HSV2 = rgb2hsv(RGB2);
            V12 = .5*(HSV1(:,:,3)+HSV2(:,:,3));
            [HSVX, HSVY, HSVZ] = sph2cart(RGBhsv(:,:,1)*2*pi,(pi/2)-sscale*RGBhsv(:,:,2),V12);
            
            % figure,
            % quiver3(zeros(100),zeros(100),zeros(100),HSVX,HSVY,HSVZ);
            
            % color displayed
            mapper = rgb2hsv(RGB/2);
            mapper2 = mapper;
            mapper2(:,:,3) = mapper(:,:,3)*2;
            RGBmapper = hsv2rgb(mapper2);
            
            mapper = hsv2rgb(cat(3,RGBhsv(:,:,1:2),V12));
            
            surf(HSVX,HSVY,HSVZ,mapper); shading flat
            xlim([-1 1]);
            ylim([-1 1]);
            hold on
            axis off
            view([a e])
            pause(0.01)
            a=a+5;
            e=40*sin(i/15) +45;
            ax=gca;
            ax.DataAspectRatio = [2 2 1];
            axis vis3d;
        end
    end % for HSV space
    
else if strcmp(play,'xxy')
        %% physical colours (wavelengths) in xyY space with dp2xy (380-780)
        
        windex = [450:5:680];
        %windex([2 4 36 38 40:end-1])=[];
        
        waves = windex;
        waves1 = [waves waves waves];
        waves2 = [circshift(waves,20) circshift(waves,30) circshift(waves,2)];
        
        a = 130;
        e = 0;
        v=1;
        figure('Color',[0 0 0]);
        for i = 1:length(waves1)
            Ys = linspace(0,v,100); % varying luminence vector
            Y2 = meshgrid(Ys,Ys);
            
            H1 = dp2xy([waves1(i) .7]); % put nms into xy chrom
            H2 = dp2xy([waves2(i) .7]); % put nms into xy chrom ,ones(1,3)
            
            h1x = repmat(H1(1),100);
            h1y = repmat(H1(2),100);
            h1scale = repmat(Ys/v,[100 1]);
            
            h2x = repmat(H2(1),100);
            h2y = repmat(H2(2),100);
            h2scale = repmat((Ys/v)',[1 100]);
            
            [h1X,h1Y,h1Z]=xy2xyz(h1x,h1y); % into tristimulus XYZ
            [h2X,h2Y,h2Z]=xy2xyz(h2x,h2y); % into tristimulus XYZ
            
            X = (h1X.*h1scale)+(h2X.*h2scale); % combine stimulus values
            Y = (h1Y.*h1scale)+(h2Y.*h2scale);
            Z = (h1Z.*h1scale)+(h2Z.*h2scale);
            
            RGB = xyz2rgb(v*X,v*Y,v*Z); % imshow(RGB);
            
            % convert to chromaticity
            [x,y,Yc] = xyz2xyy(X,Y,Z);
            
            % figure, surf(X,Y,Z,RGB); shading flat
            surf(x,y,log(Y),RGB); shading flat, hold on
            xlim([0 .8])
            ylim([0 .8])
            axis off
            view([a e])
            pause(0.1)
            a=a+5;
            e=8*sin(i/15) +52;
            ax=gca;
            %         ax.DataAspectRatio = [2 2 1];
            %         if i*j == 1
            %             camzoom(1)
            %         end
            axis vis3d;
        end
        
    elseif strcmp(play,'saturated')
        %% view all saturated colors
        figure('Color',[0 0 0]);
        HH1 = linspace(0,.99,100); % across all hues 0-1
        HH2 = HH1+.005;
        HH2 = circshift(HH1,1);
        sats = 1;
        for j = 1:length(sats)
            for i = 1:length(HH1)
                % define colors in HSV
                H1 = hsv2rgb([HH1(i),sats(j),1]); % red
                H2 = hsv2rgb([HH2(i),sats(j),1]); % yellow
                
                h1r = linspace(0,H1(1),100);
                h1g = linspace(0,H1(2),100);
                h1b = linspace(0,H1(3),100);
                h2r = linspace(0,H2(1),100);
                h2g = linspace(0,H2(2),100);
                h2b = linspace(0,H2(3),100);
                
                [RX,RY] = meshgrid(h1r,h2r);
                [GX,GY] = meshgrid(h1g,h2g);
                [BX,BY] = meshgrid(h1b,h2b);
                
                R = (RX + RY); % add components of colors
                G = (GX + GY);
                B = (BX + BY);
                
                RGB = cat(3,R,G,B);
                RGBhsv = rgb2hsv(RGB);
                [HSVX, HSVY, HSVZ] = sph2cart(RGBhsv(:,:,1)*2*pi,(pi/2)-RGBhsv(:,:,2),0.5*RGBhsv(:,:,3));
                
                % need to correct V
                
                RGB1 = cat(3,RX,GX,BX); RGB2 = cat(3,RY,GY,BY);
                HSV1 = rgb2hsv(RGB1); HSV2 = rgb2hsv(RGB2);
                V12 = .5*(HSV1(:,:,3)+HSV2(:,:,3));
                [HSVX, HSVY, HSVZ] = sph2cart(RGBhsv(:,:,1)*2*pi,(pi/2)-RGBhsv(:,:,2),V12);
                
                % color displayed
                mapper = rgb2hsv(RGB/2);
                mapper2 = mapper;
                mapper2(:,:,3) = mapper(:,:,3)*2;
                RGBmapper = hsv2rgb(mapper2);
                
                mapper = hsv2rgb(cat(3,RGBhsv(:,:,1:2),V12));
                
                surf(HSVX,HSVY,HSVZ,RGB); shading flat
                hold on
                axis off
                view([v 0])
                pause(0.1)
                v=v+10;
                ax=gca;
                ax.DataAspectRatio = [3 3 1];
                xlim([-1 1]);
                ylim([-1 1]);
                axis vis3d;
            end
        end
    elseif strcmp(play,'vectors')
        %% just test combining 2 colors vectors..
        
        sam = [0 1 0];
        mike = [0.5 0 .5];
        smike = (sam+mike)/2;
        
        [Hsam,Ssam,Vsam] = rgb2hsv(sam);
        [Hmike,Smike,Vmike] = rgb2hsv(mike);
        [Hsmike,Ssmike,Vsmike] = rgb2hsv(smike);
        
        %Vsam = Vsam/2; Vmike = Vmike/2;
        Vsmike = (Vsam+Vmike)/2; % value has to be mean
        
        [Xsam,Ysam,Zsam] = sph2cart(Hsam,Ssam,Vsam);
        [Xmike,Ymike,Zmike] = sph2cart(Hmike,Smike,Vmike);
        [Xsmike,Ysmike,Zsmike] = sph2cart(Hsmike*2*pi,(pi/2)-Ssmike,Vsmike);
        
        % figure,
        % quiver3([0 0 0],[0 0 0],[0 0 0],[Xsam Xmike Xsmike],[Ysam Ymike Ysmike],[Zsam Zmike Zsmike]);
        
        %figure
        quiver3([0],[0],[0],[Xsam],[Ysam],[Zsam]); hold on
        quiver3([0],[0],[0],[Xmike],[Ymike],[Zmike]);
        quiver3([0],[0],[0],[Xsmike],[Ysmike],[Zsmike]);
    end
end
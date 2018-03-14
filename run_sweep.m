 

clear all   % clear variables
close all   % close all figures
clc         % clear command window



%% PATH stuff for Mac, to be adapted or commented

% gmsh and getdp executables (with correct path)
%gmsh_exe  = '/Users/johangyselinck/src/gmsh/bin/gmsh' ;
%getdp_exe = '/Users/johangyselinck/src/getdp/bin/getdp' ;
%putty = 'C:\Program Files (x86)\PuTTY\putty.exe -ssh niyomsatian@192.168.0.57:niyomsatian';
%gmsh_exe  = 'C:\gmsh_getdp\gmsh-getdp-Windows64\gmsh'  ;
%getdp_exe = 'C:\gmsh_getdp\gmsh-getdp-Windows64\getdp'  ;
gmsh_exe  = 'D:\gmsh_getdp\onelab-Windows64\gmsh'  ;
getdp_exe = 'D:\gmsh_getdp\onelab-Windows64\getdp'  ;

% directory of input files
%%% file_path = '' ; 

% directory in which postprocessing files are written
Dir_res = 'coeff/' ; 


%% PATH stuff for Windows, possibly to be uncommented and adapted 

%getdp_exe = 'C:\gmsh_getdp\gmsh-getdp-Windows64\getdp ' ;
%gmsh_exe  = 'C:\gmsh_getdp\gmsh-getdp-Windows64\gmsh ' ;
%dir_res = 'res/' ; 
%%% file_path =' C:\Users\Lil-Chris\Desktop\gmsh\one_conductor_round_24feb2014\'; 





%% arguments for gmsh and getdp for the problem at hand
%msh_arg  = ' TR_ETD59_2D.geo -2' ; % make a 2D mesh (without opening the graphical interface)
%gmsh_arg  = ' cell.geo '; 
getdp_arg = ' tfo3d.pro -msh tfo3d_builtin.msh -sol Analysis -pos Get_GlobalQuantities ' ;



%% adjust RMD via 'RMD.geo', uncomment the Include in the mean geo file if need be!!
% % if 0
    % % RMD = 6 ;
    % % fid = fopen('RMD.geo','w');
    % % fprintf(fid, 'RMD = %f ; \n', RMD);
    % % fclose(fid);
% % end


% %% solving the problem
% if 1
    % % run gmsh and getdp
    % system([gmsh_exe gmsh_arg])
    % system([getdp_exe getdp_arg])
% end


% %% b vector along a horizontal line and along a circle
% % and verification of Ampere's law

% % b vector along a horizontal line
% if 0    
    % %load res/b_line.dat ; 
    % b_line = load([Dir_res 'b_line.dat']) ;
    % x_v = b_line(:,1) ;
    % y_v = b_line(:,2) ;    
    % bx_v = b_line(:,4) ;
    % by_v = b_line(:,5) ;
        
    % Xfac = 1e3 ; % mm instead of m
    % Bfac = 1e3 ; % mT insteda of T
    % figure    
    % plot(Xfac*x_v,bx_v*Bfac,'r')    
    % hold on 
    % xlabel('x (mm)')
    % grid on
    % ylabel('flux density (mT)')
    % plot(Xfac*x_v,by_v*Bfac,'b')
    % legend(' bx',' by')
% end


% % b vector along a circle
% if 0    
    % b_circle = load([Dir_res 'b_circle.dat']) ;  
    % x_v = b_circle(:,1) ;
    % y_v = b_circle(:,2) ;    
    % bx_v = b_circle(:,4) ;
    % by_v = b_circle(:,5) ;
    
    % % to be completed
    % very_wrong_v = rand(size(x_v)) ;
    % theta_v = very_wrong_v ;
    % r_v = very_wrong_v ;
    % br_v = very_wrong_v ;
    % bt_v = very_wrong_v ;
    
    % Tfac = 180/pi ; % degrees instead of radians
    % Bfac = 1e3 ; % mT instead of T
    % figure    
    % plot(theta_v*Tfac,bt_v*Bfac,'r') ;    
    % hold on
    % plot(theta_v*Tfac,br_v*Bfac,'b') ;
    % xlabel('angle (°)')
    % grid on
    % ylabel('flux density (mT)')
    % legend(' tangential component',' radial component', 'Location','Best')
    % title('Flux density on a circle')
   
    % % verificaton of Ampere's law (satisfied in a weak, i.e. averaged, way)     
    % mu0 = 4*pi*1e-7 ;
    % MMF = mean(bt_v)* r_v(1)*2*pi/mu0         
% end


%% energy versus RMD, for checking convergence

if 1
    %for a = [0.00019 0.0002 0.00021] % 0.1e-3*2*(0.95 ... 1.05)
        %for b = [0.00019 0.0002 0.00021]%[0.0002 0.00021]%[0.00019 0.0002 0.00021]
    
    %RMD_v = 1:0.5:5 ; % vector with relative mesh density (RMD) values
    
    %RMD_v = sqrt(2).^(0:8) ; % exponentially increasing
    
    %RMD_v = 5
    
	%fill_v = 0.5:0.5:2 ; % vector with fill factor values 
    Rc = 0.1e-3;%0.05;
	X_v = linspace(0.2,2,10); % vector with frequency values 
	
	hex = 1; square = 0;
	skin = 1; prox = 2;
    
    r = 0.1e-3; % strand radius
    sigma = 4.35e7;
    mu0 = 4*pi*1e-7;
    
    Irms_pri = 1;  ph_pri = 0;      % 300 150 0
    Irms_sec = 0;    ph_sec = 0;  % 300 100 100 180
    %strand_dia_pri = a;
    %strand_dia_sec = b;
    %         set RMD = 1.5
    
%     dum = load([Dir_res 'pB_RS_la0.38.dat']) ;
%     pB_pri = dum(:,2);
%     dum = load([Dir_res 'qB_RS_la0.38.dat']) ;
%     qB_pri = dum(:,2);
%     dum = load([Dir_res 'pI_RS_la0.38.dat']) ;
%     pI_pri = dum(:,2);
%     dum = load([Dir_res 'qI_RS_la0.38.dat']) ;
%     qI_pri = dum(:,2);
%     
%     dum = load([Dir_res 'pB_RS_la0.63.dat']) ;
%     pB_sec = dum(:,2);
%     dum = load([Dir_res 'qB_RS_la0.63.dat']) ;
%     qB_sec = dum(:,2);
%     dum = load([Dir_res 'pI_RS_la0.63.dat']) ;
%     pI_sec = dum(:,2);
%     dum = load([Dir_res 'qI_RS_la0.63.dat']) ;
%     qI_sec = dum(:,2);

%     fid = fopen('param.pro','w');
%     fprintf(fid, 'pB_pri_00 = %f ; \n', pB_pri);
%     fprintf(fid, 'qB_pri_00 = %f ; \n', qB_pri);
%     fprintf(fid, 'pI_pri_00 = %f ; \n', pI_pri); 
%     fprintf(fid, 'qI_pri_00 = %f ; \n', qI_pri);
%     
%     fprintf(fid, 'pB_sec_00 = %f ; \n', pB_sec);
%     fprintf(fid, 'qB_sec_00 = %f ; \n', qB_sec);
%     fprintf(fid, 'pI_sec_00 = %f ; \n', pI_sec); 
%     fprintf(fid, 'qI_sec_00 = %f ; \n', qI_sec);
%     
%     fprintf(fid, 'Ic_pk_pri_00 = %f ; \n', Ipk_pri);
%     fprintf(fid, 'Ic_ph_pri_00 = %f ; \n', ph_pri);
%     fprintf(fid, 'Ic_pk_sec_00 = %f ; \n', Ipk_sec);
%     fprintf(fid, 'Ic_ph_sec_00 = %f ; \n', ph_sec);
%     
%     %fclose(fid);    
%     system([gmsh_exe gmsh_arg])
    
        f_v = logspace(1,6,20);
        for ix = 1:length(f_v)%1:length(X_v)
            tic        
            %X = X_v(ix)
            %X = 2;
            %X = 0.1

            %f = 100e3;

            %f = (X/r)^2/(pi*sigma*mu0);
            f = f_v(ix);


            % adjust RMD via 'RMD.geo'
            %fid = fopen('cell_param.pro','w');
            %fprintf(fid, 'RMD = %f ; \n', RMD);
            %fid = fopen('cell_param.pro','w');

            fid = fopen('param_init.pro','w');
            fprintf(fid, 'Freq_00 = %f ; \n', f);
    %         fprintf(fid, 'pB_pri_00 = %f ; \n', pB_pri(ix));
    %         fprintf(fid, 'qB_pri_00 = %f ; \n', qB_pri(ix));
    %         fprintf(fid, 'pI_pri_00 = %f ; \n', pI_pri(ix)); 
    %         fprintf(fid, 'qI_pri_00 = %f ; \n', qI_pri(ix));
    % 
    %         fprintf(fid, 'pB_sec_00 = %f ; \n', pB_sec(ix));
    %         fprintf(fid, 'qB_sec_00 = %f ; \n', qB_sec(ix));
    %         fprintf(fid, 'pI_sec_00 = %f ; \n', pI_sec(ix)); 
    %         fprintf(fid, 'qI_sec_00 = %f ; \n', qI_sec(ix));

            fprintf(fid, 'Irms_pri_00 = %f ; \n', Irms_pri);
            %fprintf(fid, 'Ic_ph_pri_00 = %f ; \n', ph_pri);
            fprintf(fid, 'Irms_sec_00 = %f ; \n', Irms_sec);
            %fprintf(fid, 'Ic_ph_sec_00 = %f ; \n', ph_sec);
            %fprintf(fid, 'strand_dia_pri_00 = %f ; \n', strand_dia_pri);
            %fprintf(fid, 'strand_dia_sec_00 = %f ; \n', strand_dia_sec);

            fclose(fid);

            % run gmsh and getdp
            %system([gmsh_exe gmsh_arg])
            %%% run gmsh seperately at the finest mesh !!!!!!!
            system([getdp_exe getdp_arg])

             % read output files
    % 		dum = load([Dir_res 'pB_RS_la0.2.dat']) ;
    %         Xp(ix) = dum(ix,1);
    % 		p(ix) = dum(ix,2);
    %         
    %         dum_2 = load([Dir_res 'qB_RS_la0.2.dat']) ;
    %         Xq(ix) = dum_2(ix,1);
    %         q(ix) = dum_2(ix,2);
    %         %dum = load([Dir_res 'MagEnergy.dat']) ;
    %         %MagEnergy_v(irmd) = dum(2)

            calc_time_v(ix) = toc % calculation time per RMD value       
        end
       %end
    %end
  
end


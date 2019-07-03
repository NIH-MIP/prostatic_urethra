%% DESCRIPTION
%
% This program takes in T2 MRI data along with .voi segmentation files of 
% urethra, TZ, and PZ zones. It decimates the urethra segmentation to a 
% center of mass volume, from which a piecewise spline is fit to the data
% and the track-length as well as urethra angles are calculated.
% please note multiple methods for angle calculation are used in this code
% but only one was used for the results in publication by Sanford et al
% which correpnds to "infl_angle_dist" in the output sturctures
%
% 3D mesh created for visualization purposes
%
% inputs:
%    - voiLoc : urethra segmentation location/filename
%    - t2Loc : directory where t2 files are located
%    - tzLoc : tz segmentation location/filename
%    - wpLoc : wp segmentation location/filename
%    - saveDir : directory for output file (in the form of matlab structure file
%
% outputs:
%   - matlab structure file containing all data 
%             - t2 folder
%             - voi_id
%             - wp_id
%             - tz_id 
%             - voxel size (of t2 dicom)
%             - t2_data : structure contianing full stack of dicom images
%             - voi_points : urethra points in each T2 slice
%             - com_points : center of mass in each T2 slice
%             - com_shp : alpha shape of mesh
%             - poly : polygon of urethra com-volume points
%             - sag: structure containing
%                     - infl_ind : index of inflection point from derivatives
%                     - infl_ind_dist : index of inflection point from deviation (USED IN MANUSCRIPT)
%                     - length base2apex straight (straight line)
%                     - length base2space piecewise (piece along splinefit points)
%                     - infl_angle : angle of inflection point from derivatives
%                     - infl_angle_dist : angle of inflection point from deviation (USED IN MANUSCRIPT)
%                     - tz_prop : proportion of tz volume above plane from derivatives
%                     - tz_prop : proportion of tz volume above plane from deviation (USED IN MANUSCRIPT)
%                     - wp_prop : proportion of wp volume above plane from derivatives
%                     - wp_prop : proportion of wp volume above plane from deviation (USED IN MANUSCRIPT)
%             - cor: see above
%             - tz volume
%             - wp volume

%   - there is an area in the code that can be uncommented to output as csv too       
%% DEPENDENCIES
%
% - local installation of MATLAB 
%
% - MATLAB Image Processing Toolbox
%
% - spline fit function from 
%     https://www.mathworks.com/matlabcentral/fileexchange/71225-splinefit   
%

%%
function [] = Urethra_Characterization(voiLoc,t2Loc,tzLoc,wpLoc,saveDir)
%% read in voi and t2
    voi_points = VOIreader_fp(voiLoc); %re-wrote to accept floating point values for voxel point
    wp_points = VOIreader_fp(wpLoc);
    if(wp_points{1,1} == 0)
        wp_points = wp_points(2:end,:);
    end
    tz_points = VOIreader_fp(tzLoc);
    if(tz_points{1,1} == 0)
        tz_points = tz_points(2:end,:);
    end
    t2 = ConvertDicom(t2Loc); %grab image resolution info from t2
    t2parts = strsplit(t2Loc,filesep);
    pt_id = t2parts{find(strcmpi(t2parts,'dicoms'))-1};
    mkdir([saveDir filesep pt_id])

% find points of voi: floating point and mask

    all_pts = []; %this is a list of all polygon points contained in VOI
    mask = zeros(size(t2.data)); %this will be full mask from VOI
    mask_sm = zeros(size(t2.data)); %this will be COM mask with volumetric spread
    com_pts = zeros(size(voi_points,1),3); %this will be exact COM points

    %vois are built on each slice
    for i = 1:size(voi_points,1)
        slice_i = voi_points{i,1}; %slice level
        pts_i = voi_points{i,2};   %polygon points from voi
        slice_mask = poly2mask(pts_i(:,1),pts_i(:,2),size(t2.data,1),size(t2.data,2)); %create mask from points
        mask(:,:,slice_i) = slice_mask; %fill mask at every slice

        %find centroid (COM) of polygon on slice level
        s = regionprops(slice_mask,slice_mask,'Centroid');
        com_pts(i,1) = s.Centroid(1); 
        com_pts(i,2) = s.Centroid(2);
        com_pts(i,3) = slice_i;

        %fit a circular mask to COM points where radius = [resolution z plane]/[reslution x-y plane]
        radii = (t2.voxel_size(3)/t2.voxel_size(1))/2;
        xDim = size(slice_mask,1);
        yDim = size(slice_mask,2);
        [xx,yy] = meshgrid(1:yDim,1:xDim); %remember plotting vs matrices are different x-y indexing
        mask_shrink = false(xDim,yDim);
        mask_shrink = mask_shrink | hypot(xx - round(s.Centroid(1)), yy - round(s.Centroid(2))) <= radii;
        mask_sm(:,:,slice_i) = mask_shrink;

        pts_i(:,3) = slice_i; %keep everything in voxel coordinates for now 
        all_pts = cat(1, all_pts, pts_i);
    end

    %wp mask
    wp_mask = zeros(size(t2.data)); %this will be full mask from VOI
    wp_pts = [];
    for i = 1:size(wp_points,1)
        slice_i = wp_points{i,1}; %slice level
        pts_i = wp_points{i,2};   %polygon points from voi
        slice_mask = poly2mask(pts_i(:,1),pts_i(:,2),size(t2.data,1),size(t2.data,2)); %create mask from points
        wp_mask(:,:,slice_i) = slice_mask; %fill mask at every slice
        pts_i(:,3) = slice_i;
        wp_pts = cat(1, wp_pts, pts_i);
    end
    
    %tz mask
    tz_mask = zeros(size(t2.data)); %this will be full mask from VOI
    tz_pts = [];
    for i = 1:size(tz_points,1)
        slice_i = tz_points{i,1}; %slice level
        pts_i = tz_points{i,2};   %polygon points from voi
        slice_mask = poly2mask(pts_i(:,1),pts_i(:,2),size(t2.data,1),size(t2.data,2)); %create mask from points
        tz_mask(:,:,slice_i) = slice_mask; %fill mask at every slice
        pts_i(:,3) = slice_i;
        tz_pts = cat(1, tz_pts, pts_i);
    end    
    
% 3D renderings in real space

    %first grab all polygon points from voi file and multiply by spatial resolution 
    P_full(:,1) = all_pts(:,2)*t2.voxel_size(1);
    P_full(:,2) = all_pts(:,1)*t2.voxel_size(2);
    P_full(:,3) = all_pts(:,3)*t2.voxel_size(3);
    %fit a simple triangulation with shrink factor 1 (tightest possible combinations)
    %note to self: we could do better than this, but this step is mainly
    %for visualization so not important to consider at the moment
    k = boundary(P_full,1); 
    %now plot all points from this boundary object 
    figure, plot3(P_full(:,1),P_full(:,2),P_full(:,3),'.','MarkerSize',10)
    trisurf(k,P_full(:,1),P_full(:,2),P_full(:,3),'Facecolor','blue','FaceAlpha',0.01)
    full_polygon = struct;
    full_polygon.P_full = P_full;
    full_polygon.k = k;

    %next grab all points contained within COM polygon mask
    findpts = find(mask_sm>0);
    [Px,Py,Pz] = ind2sub(size(t2.data),findpts); 
    P = cat(2,Px, Py, Pz);
    hold on
    %create alpha shape from COM polygon mask using the SAME radius we used to generate the circular masks. alpha shapes create convex hull
    %changes from a set of points (sub of delaunay triangulation) using the same alpha radius is important here as it should connect to
    %one volumetric region following the COM volumetric along the entire prostate
    shp = alphaShape([P(:,1)*t2.voxel_size(1),P(:,2)*t2.voxel_size(2),P(:,3)*t2.voxel_size(3)],(t2.voxel_size(3)/t2.voxel_size(1)));
    plot(shp,'Facecolor','red','FaceAlpha',0.3)
    title('Surface mesh of original VOI and COM-region')

    
    P_wp(:,1) = wp_pts(:,2)*t2.voxel_size(1);
    P_wp(:,2) = wp_pts(:,1)*t2.voxel_size(2);
    P_wp(:,3) = wp_pts(:,3)*t2.voxel_size(3);
    k_wp = boundary(P_wp,1);
    %trisurf(k_wp,P_wp(:,1),P_wp(:,2),P_wp(:,3),'Facecolor','yellow','FaceAlpha',0.05,'EdgeAlpha','0.1')
    
    P_tz(:,1) = tz_pts(:,2)*t2.voxel_size(1);
    P_tz(:,2) = tz_pts(:,1)*t2.voxel_size(2);
    P_tz(:,3) = tz_pts(:,3)*t2.voxel_size(3);
    k_z = boundary(P_tz,1);
    %trisurf(k_tz,P_tz(:,1),P_tz(:,2),P_tz(:,3),'Facecolor','green','FaceAlpha',0.05,'EdgeAlpha','0.1')
   
    % find piece-wise linear distance along all COMS
    pw_length = 0;
    for j = 2:size(com_pts,1)
        v1 = [com_pts(j-1,1)*t2.voxel_size(1), com_pts(j-1,2)*t2.voxel_size(2), com_pts(j-1,3)*t2.voxel_size(3)];
        v2 = [com_pts(j,1)*t2.voxel_size(1), com_pts(j,2)*t2.voxel_size(2), com_pts(j,3)*t2.voxel_size(3)];
        length_j = norm(v1-v2);
        pw_length = pw_length + length_j;
    end
    length3D_base2apex_piecewise = pw_length;
    vstart = [com_pts(1,1)*t2.voxel_size(1), com_pts(1,2)*t2.voxel_size(2), com_pts(1,3)*t2.voxel_size(3)];
    vend = [com_pts(end,1)*t2.voxel_size(1), com_pts(end,2)*t2.voxel_size(2), com_pts(end,3)*t2.voxel_size(3)];
    length3D_base2apex_straight = norm(vstart-vend);
    
%% summary stats
    ustats = struct;
    ustats.t2_folder = t2Loc;
    ustats.voi_id = voiLoc;
    ustats.wp_id = wpLoc;
    ustats.tz_id = tzLoc;
    ustats.voxel_size = t2.voxel_size;
    ustats.t2_data = t2;
    ustats.voi_points = voi_points;
    ustats.com_points = com_pts;
    ustats.com_shp = shp;
    ustats.poly = full_polygon;
    ustats.length3D_base2apex_piecewise = length3D_base2apex_piecewise;
    ustats.length3D_base2apex_straight = length3D_base2apex_straight;
    
    sag = inflection_calculator(P, com_pts,t2, 'sagittal');
     view(90,-90)
    
    cor = inflection_calculator(P, com_pts,t2, 'coronal');
     view(90,-90)
    
    ustats.sag = sag;
    ustats.cor = cor;
    
    tz_inds = find(tz_mask>0);
    wp_inds = find(wp_mask>0);
    
    tz_vol = length(tz_inds)*t2.voxel_size(1)*t2.voxel_size(2)*t2.voxel_size(3)*1e-3;
    wp_vol = length(wp_inds)*t2.voxel_size(1)*t2.voxel_size(2)*t2.voxel_size(3)*1e-3;
    
    ustats.tz_vol = tz_vol;
    ustats.wp_vol = wp_vol;
    
    %proportion of tz and wp for inflection angle
    if(ustats.sag.infl_angle < 180)
        coefficients = polyfit([ustats.sag.infl_ind(3), ustats.com_points(end,3)],[ustats.sag.infl_ind(2), ustats.com_points(end,2)] , 1);
        plane_mask = zeros(size(t2.data));
        for i= 1:size(t2.data,3)
            z_val = i;
            s_val = round(z_val*coefficients(1)+coefficients(2));
            if(1 <= s_val && s_val <= size(t2.data,1))
                plane_mask(:,s_val:end,z_val) = 1;
            end
        end
        plane_inds = find(plane_mask>0);
        tz_prop = length(intersect(tz_inds,plane_inds))/length(tz_inds);
        wp_prop = length(intersect(wp_inds,plane_inds))/length(wp_inds);
        ustats.sag.tz_prop = tz_prop;
        ustats.sag.wp_prop = wp_prop;
    end
    
    %proportion of tz and wp for max dist angle
        coefficients = polyfit([ustats.sag.infl_ind_dist(3), ustats.com_points(end,3)],[ustats.sag.infl_ind_dist(2), ustats.com_points(end,2)] , 1);
        plane_mask = zeros(size(t2.data));
        for i= 1:size(t2.data,3)
            z_val = i;
            s_val = round(z_val*coefficients(1)+coefficients(2));
            if(1 <= s_val && s_val <= size(t2.data,1))
                plane_mask(:,s_val:end,z_val) = 1;
            end
        end
        plane_inds = find(plane_mask>0);
        tz_prop = length(intersect(tz_inds,plane_inds))/length(tz_inds);
        wp_prop = length(intersect(wp_inds,plane_inds))/length(wp_inds);
        ustats.sag.tz_prop_dist = tz_prop;
        ustats.sag.wp_prop_dist = wp_prop;

% if you want to save everything out to text file, uncomment below
%     save([saveDir filesep pt_id filesep pt_id],'ustats','-v7.3');
%     
%     fileID = fopen([saveDir filesep pt_id filesep pt_id '.csv'],'w');
%     
%     %volume stats
%     fprintf(fileID,'%s, %s \n','tz_vol',num2str(ustats.tz_vol));
%     fprintf(fileID,'%s, %s \n','wp_vol',num2str(ustats.wp_vol));
%     
%     %length stats
%     fprintf(fileID,'%s, %s \n','sag_length_st',num2str(ustats.sag.length_base2apex_straight));
%     fprintf(fileID,'%s, %s \n','sag_length_pc',num2str(ustats.sag.length_base2apex_piecewise));
%     fprintf(fileID,'%s, %s \n','sag_length_ratio',num2str(ustats.sag.length_base2apex_piecewise/ustats.sag.length_base2apex_straight));
%     fprintf(fileID,'%s, %s \n','cor_length_st',num2str(ustats.cor.length_base2apex_straight));
%     fprintf(fileID,'%s, %s \n','cor_length_pc',num2str(ustats.cor.length_base2apex_piecewise));
%     fprintf(fileID,'%s, %s \n','cor_length_ratio',num2str(ustats.cor.length_base2apex_piecewise/ustats.cor.length_base2apex_straight));
%     fprintf(fileID,'%s, %s \n','3d_length_st',num2str(ustats.length3D_base2apex_straight));
%     fprintf(fileID,'%s, %s \n','3d_length_pc',num2str(ustats.length3D_base2apex_piecewise));
%     fprintf(fileID,'%s, %s \n','3d_length_ratio',num2str(ustats.length3D_base2apex_piecewise/ustats.length3D_base2apex_straight));
%     
%     %angles (inflection and max distance)
%     fprintf(fileID,'%s, %s \n','sag_angle_infl',num2str(ustats.sag.infl_angle));
%     fprintf(fileID,'%s, %s \n','cor_angle_infl',num2str(ustats.cor.infl_angle));
%     fprintf(fileID,'%s, %s \n','sag_angle_maxdist',num2str(ustats.sag.infl_angle_dist));
%     fprintf(fileID,'%s, %s \n','cor_angle_maxdist',num2str(ustats.cor.infl_angle_dist));
%     
%     %locations (inflection and max distance)
%     fprintf(fileID,'%s, %s \n','urethra_apex_slice',num2str(ustats.com_points(1,3)));
%     fprintf(fileID,'%s, %s \n','urethra_base_slice',num2str(ustats.com_points(end,3)));
%     fprintf(fileID,'%s, %s \n','urethra_sag_slice_infl',num2str(ustats.sag.infl_ind(3)));
%     fprintf(fileID,'%s, %s \n','urethra_cor_slice_infl',num2str(ustats.cor.infl_ind(3)));
%     fprintf(fileID,'%s, %s \n','urethra_sag_slice_dist',num2str(ustats.sag.infl_ind_dist(3)));
%     fprintf(fileID,'%s, %s \n','urethra_cor_slice_dist',num2str(ustats.cor.infl_ind_dist(3)));
%     
%     if(ustats.sag.infl_angle < 180)
%         fprintf(fileID,'%s, %s \n','infl_sag_tz_ratio',num2str(ustats.sag.tz_prop));
%         fprintf(fileID,'%s, %s \n','infl_sag_tz_vol_above',num2str(ustats.sag.tz_prop*ustats.tz_vol));
%         fprintf(fileID,'%s, %s \n','infl_sag_tz_vol_below',num2str((1-ustats.sag.tz_prop)*ustats.tz_vol));
%         fprintf(fileID,'%s, %s \n','infl_sag_wp_ratio',num2str(ustats.sag.wp_prop));
%         fprintf(fileID,'%s, %s \n','infl_sag_wp_vol_above',num2str(ustats.sag.wp_prop*ustats.wp_vol));
%         fprintf(fileID,'%s, %s \n','infl_sag_wp_vol_below',num2str((1-ustats.sag.wp_prop)*ustats.wp_vol));
%     else
%         fprintf(fileID,'%s, %s \n','infl_sag_tz_ratio','NA');
%         fprintf(fileID,'%s, %s \n','infl_sag_tz_vol_above','NA');
%         fprintf(fileID,'%s, %s \n','infl_sag_tz_vol_below','NA');
%         fprintf(fileID,'%s, %s \n','infl_sag_wp_ratio','NA');
%         fprintf(fileID,'%s, %s \n','infl_sag_wp_vol_above','NA');
%         fprintf(fileID,'%s, %s \n','infl_sag_wp_vol_below','NA');
%     end
%     
%     fprintf(fileID,'%s, %s \n','dist_sag_tz_ratio',num2str(ustats.sag.tz_prop_dist));
%     fprintf(fileID,'%s, %s \n','dist_sag_tz_vol_above',num2str(ustats.sag.tz_prop_dist*ustats.tz_vol));
%     fprintf(fileID,'%s, %s \n','dist_sag_tz_vol_below',num2str((1-ustats.sag.tz_prop_dist)*ustats.tz_vol));
%     fprintf(fileID,'%s, %s \n','dist_sag_wp_ratio',num2str(ustats.sag.wp_prop_dist));
%     fprintf(fileID,'%s, %s \n','dist_sag_wp_vol_above',num2str(ustats.sag.wp_prop_dist*ustats.wp_vol));
%     fprintf(fileID,'%s, %s \n','dist_sag_wp_vol_below',num2str((1-ustats.sag.wp_prop_dist)*ustats.wp_vol));
%     
%     
%     fclose(fileID);
end


%%
function planestats = inflection_calculator(P, com_pts, t2, plane)
    planestats = struct;
    planestats.plane = plane;
    switch plane
        case 'sagittal'
            com_column = 2;
            P_column  = 1;
        case 'coronal'
            com_column = 1;
            P_column  = 2;
    end

    y = P(:,P_column)*t2.voxel_size(P_column); %remember P contains the indices of all points in COM polygon
    x = P(:,3)*t2.voxel_size(3); %note we are now in resolution-adjusted space, so we can measure distances
    
    %we will set constraints on our piecewise spline function to begin and
    %end at COM points from the first and last slice of urethra segmentations
    yc = [round(com_pts(1,com_column)*t2.voxel_size(com_column)),round(com_pts(end,com_column)*t2.voxel_size(com_column))]; %first and last centroid x points
    xc = [round(com_pts(1,3)*t2.voxel_size(3)),round(com_pts(end,3)*t2.voxel_size(3))]; %first and last centroid z points
    con = struct('xc',xc,'yc',yc);
    
    %breaks are points that specify the piecewise segments, exact COM points in x points is a natural choice 
    %breaks = linspace(xc(1),xc(end),round((cc(end)-cc(1))/t2.voxel_size(3)));
    breaks = [round(com_pts(:,3)*t2.voxel_size(3))];
    
    %first plot com polygon points 
    figure, plot(x,y,'.')
    title([plane ' plane: spline fit from COM-region points'])
    xx = linspace(xc(1),xc(end),400); %400 x points
    
    %fit spline
    pp1 = splinefit(x,y,3,con);
    
    y1 = ppval(pp1,xx); % sampled 400 z points (NOTE WE MIGHT NOT HAVE GREAT FIT, DONT WORRY YET)
    qq = ppdiff(pp1,2); % now take second order derivative to identify inflectin points
    qq_vals = ppval(qq,xx); %grab those 2nd order vals

    %here we are interested in two things:
    %1) minimum of the 2nd derivative to find most concave down region
        infl_ind = find(qq_vals == min(qq_vals));
    %2) any change from concave to convex (i.e. wiggles on the plot)    
        infl_change = find(diff(qq_vals>=0,1)~=0); %yes we did just take derivative of 2nd order derivative values

    %find x points associated with the inflection point and wiggle points    
    x_infl = xx(infl_ind);
    x_change = xx(infl_change);
    %also grab y values (NOTE THESE ARE THE SAMPLED DATA)
    y_infl = y1(infl_ind);
    y_change = y1(infl_change);
    %and 2nd-O Deriv values for funsies
    qq_infl = qq_vals(infl_ind);
    qq_change = qq_vals(infl_change);

    %plot piecewise spline fit and show inflection point + wiggle points
    hold on
    plot(xx,y1)
    plot(x_infl,y_infl,'r*')
    %plot(x_change,y_change,'b*')
     
    % fit line from this area and find orthogonal distances
    coefficients = polyfit([xc(1), xc(2)],[yc(1), yc(2)] , 1);
    
    testline = coefficients(1)*xx+coefficients(2);
    
    v1 = [xc(1), yc(1),0];
    v2 = [xc(2), yc(2),0];
    a = v1-v2;
    d_vals = zeros(size(com_pts,1),3);
    for pointi = 1:size(d_vals,1)
        p1_x = round(com_pts(pointi,3)*t2.voxel_size(3));
        p1_y = ppval(pp1,p1_x);
        p1 = [p1_x, p1_y, 0];
        
        %distance 
        l_a = norm(p1-v2);
        l_b = norm(v1-p1);
        l_c = norm(v1-v2);  
        angle_p1 = acosd((l_a*l_a+l_b*l_b-l_c*l_c)/(2*l_a*l_b));
        angle_v1 = asind((l_a*sind(angle_p1))/l_c);
        d_1 = sind(angle_v1)*l_b;
        d_2 = cosd(angle_v1)*l_b;
           
        p2_x = v1(1) + sqrt((d_2*d_2)/(1+coefficients(1)*coefficients(1)));
        p2_y = coefficients(1)*p2_x + coefficients(2);
        
        b = v2-p1;
        d = norm(cross(a,b)) / norm(a);
        d_vals(pointi,1) = d;
        d_vals(pointi,2) = p1_x;
        d_vals(pointi,3) = p1_y;
        d_vals(pointi,4) = p2_x;
        d_vals(pointi,5) = p2_y;
        plot([p1_x p2_x],[p1_y p2_y],'Color',[0.8 0.8 0.8])
        clear d p1_x p1_y p1
    end
    
    additional_infl = d_vals(find(d_vals(:,1) == max(d_vals(:,1))),2:3);
    plot(additional_infl(1),additional_infl(2),'b*')
    plot(xx,testline,'g')
    
    x_intersect = (additional_infl(1) + coefficients(1)*additional_infl(2) - coefficients(1)*coefficients(2))/(coefficients(1)*coefficients(1)+1);
    y_intersect = coefficients(1)*x_intersect + coefficients(2);
    %we have previously flipped the data, lets flip it back for finding other stats
    y_infl_true = y_infl;
    y_infl = x_infl;
    
    y_infl_true_dist = additional_infl(2);
    y_infl_dist = additional_infl(1);

    % if inflection point inside com polygon, accept it
    % if its not, find the z position its between and average between COMs 
    % report z slice of inflection point (to later map to base-mid-spex of prostate
    if((find(round(y_infl) == com_pts(:,3)*t2.voxel_size(3))>0))
        ind_pos = find(com_pts(:,3)*t2.voxel_size(3) == round(y_infl));
        ind_opts = P(find(P(:,3)*t2.voxel_size(3) == round(y_infl)),:);
        ind_val = median(ind_opts(find(round(ind_opts(:,P_column)*t2.voxel_size(1)) == round(y_infl_true)),:));
        planestats.infl_ind = [ind_val(2) ind_val(1) ind_val(3)]; %P values   
    else %if its in between levels, find the closest one
        ind_pos = find(mod(com_pts(:,3)*t2.voxel_size(3),round(y_infl)) == min(mod(com_pts(:,3)*t2.voxel_size(3),round(y_infl))));
        ind_opts = P(find(P(:,3)*t2.voxel_size(3) == com_pts(ind_pos,3)*t2.voxel_size(3)),:);
        searchvals = ind_opts(find(round(ind_opts(:,P_column)*t2.voxel_size(1)) == round(y_infl_true)),:);
        if(size(searchvals,1) == 0)
            if(round(y_infl_true) < min(round(ind_opts(:,P_column)*t2.voxel_size(1))))
                searchvals = ind_opts(find(round(ind_opts(:,P_column)*t2.voxel_size(1)) == min(round(ind_opts(:,P_column)*t2.voxel_size(1)))),:);
                ind_val = median(searchvals);
            elseif(round(y_infl_true) > max(round(ind_opts(:,P_column)*t2.voxel_size(1))))
                searchvals = ind_opts(find(round(ind_opts(:,P_column)*t2.voxel_size(1)) == max(round(ind_opts(:,P_column)*t2.voxel_size(1)))),:);
                ind_val = median(searchvals); 
            end
        else  
            ind_val = median(searchvals);
        end
        planestats.infl_ind = [ind_val(2) ind_val(1) ind_val(3)]; %P values
    end
    
    %ind_pos_dist = find(com_pts(:,3)*t2.voxel_size(3) == round(y_infl_dist));
    ind_opts = P(find(P(:,3)*t2.voxel_size(3) == round(y_infl_dist)),:);
    ind_val = median(ind_opts(find(round(ind_opts(:,P_column)*t2.voxel_size(1)) == round(y_infl_true_dist)),:));
    planestats.infl_ind_dist = [ind_val(2) ind_val(1) ind_val(3)]; %P values 
    
    % find distance between first and last centroids
    point_base = [round(com_pts(end,com_column)*t2.voxel_size(com_column)) round(com_pts(end,3)*t2.voxel_size(3))];
    point_apex = [round(com_pts(1,com_column)*t2.voxel_size(com_column)) round(com_pts(1,3)*t2.voxel_size(3))];
    length_ba = norm(point_base-point_apex);
    planestats.length_base2apex_straight = length_ba;
    
    % find piece-wise linear distance along all COMS
    pw_length = 0;
    for j = 2:size(com_pts,1)
        length_j = norm([round(com_pts(j-1,com_column)*t2.voxel_size(com_column)) round(com_pts(j-1,3)*t2.voxel_size(3))] - [round(com_pts(j,com_column)*t2.voxel_size(com_column)) round(com_pts(j,3)*t2.voxel_size(3))]);
        pw_length = pw_length + length_j;
    end
    planestats.length_base2apex_piecewise = pw_length;
    
    %fit a triangle between straight line connecting first/last centroids and inflection point
    %calculate angle between these points (hypothesis -> sharper (smaller) angle = bad)
    if(ind_pos == 1 | ind_pos == size(com_pts,1))
        angle_i  = 180;
        planestats.infl_angle = 180;
    else
        %calculate angle 
        %remember this is sagittal so we just need two points here
        point_infl = [round(planestats.infl_ind(com_column)*t2.voxel_size(com_column)) round(planestats.infl_ind(3)*t2.voxel_size(3))];
        length_bi = norm(point_base-point_infl);
        length_ai = norm(point_apex-point_infl);            
        angle_i = acosd((length_bi^2 + length_ai^2 - length_ba^2)/(2*length_bi*length_ai));
        planestats.infl_angle = angle_i;
    end

    point_infl = [round(planestats.infl_ind_dist(com_column)*t2.voxel_size(com_column)) round(planestats.infl_ind_dist(3)*t2.voxel_size(3))];
    length_bi = norm(point_base-point_infl);
    length_ai = norm(point_apex-point_infl);            
    angle_i = acosd((length_bi^2 + length_ai^2 - length_ba^2)/(2*length_bi*length_ai));
    planestats.infl_angle_dist = angle_i;
    % THIS INFLECTION ANGLE AND METHOD USED FOR PAPER
end
classdef GUI6_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        SinglefilevisualizationTab     matlab.ui.container.Tab
        TextArea_2                     matlab.ui.control.TextArea
        EditField_6                    matlab.ui.control.EditField
        FolderButton_2                 matlab.ui.control.Button
        RlowerEditField                matlab.ui.control.NumericEditField
        RlowerEditFieldLabel           matlab.ui.control.Label
        RupperEditField                matlab.ui.control.NumericEditField
        RupperEditFieldLabel           matlab.ui.control.Label
        RrotEditField                  matlab.ui.control.NumericEditField
        RrotEditFieldLabel             matlab.ui.control.Label
        DistancebetweenprotonandsamplemEditField_2  matlab.ui.control.NumericEditField
        DistancebetweenprotonandsamplemEditField_2Label  matlab.ui.control.Label
        PixelspertubeEditField_2       matlab.ui.control.NumericEditField
        PixelspertubeEditField_2Label  matlab.ui.control.Label
        NumberoftubesEditField_2       matlab.ui.control.NumericEditField
        NumberoftubesEditField_2Label  matlab.ui.control.Label
        PlotsettingsPanel              matlab.ui.container.Panel
        dmaxEditField                  matlab.ui.control.NumericEditField
        dmaxEditFieldLabel             matlab.ui.control.Label
        dminEditField                  matlab.ui.control.NumericEditField
        dminEditFieldLabel             matlab.ui.control.Label
        FileindexEditField             matlab.ui.control.NumericEditField
        FileindexEditFieldLabel        matlab.ui.control.Label
        PlotButton_2                   matlab.ui.control.Button
        UIAxes4                        matlab.ui.control.UIAxes
        PlotTab                        matlab.ui.container.Tab
        PoleFigureSettingPanel         matlab.ui.container.Panel
        GeneratePoleFigureButton       matlab.ui.control.Button
        RotatearoundYdegreesEditField  matlab.ui.control.NumericEditField
        RotatearoundYdegreesEditFieldLabel  matlab.ui.control.Label
        RotatearoundXdegreesEditField  matlab.ui.control.NumericEditField
        RotatearoundXdegreesEditFieldLabel  matlab.ui.control.Label
        PlotSettingPanel               matlab.ui.container.Panel
        ToEditField_2                  matlab.ui.control.NumericEditField
        ToEditField_2Label             matlab.ui.control.Label
        PlotdfromEditField             matlab.ui.control.NumericEditField
        PlotdfromEditFieldLabel        matlab.ui.control.Label
        NumberofpixelstoneglectattubesendEditField  matlab.ui.control.NumericEditField
        NumberofpixelstoneglectattubesendEditFieldLabel  matlab.ui.control.Label
        ResolutionDropDown             matlab.ui.control.DropDown
        ResolutionDropDownLabel        matlab.ui.control.Label
        ToEditField                    matlab.ui.control.NumericEditField
        ToEditFieldLabel               matlab.ui.control.Label
        PlotfilesfromEditField         matlab.ui.control.NumericEditField
        PlotfilesfromEditFieldLabel    matlab.ui.control.Label
        EditField                      matlab.ui.control.EditField
        FolderButton                   matlab.ui.control.Button
        PlotButton                     matlab.ui.control.Button
        ExperimentSetupPanel           matlab.ui.container.Panel
        OffsetangleEditField           matlab.ui.control.NumericEditField
        OffsetangleEditFieldLabel      matlab.ui.control.Label
        PixelspertubeEditField         matlab.ui.control.NumericEditField
        PixelspertubeEditFieldLabel    matlab.ui.control.Label
        NumberoftubesEditField         matlab.ui.control.NumericEditField
        NumberoftubesEditFieldLabel    matlab.ui.control.Label
        DistancebetweenprotonandsamplemEditField  matlab.ui.control.NumericEditField
        DistancebetweenprotonandsamplemEditFieldLabel  matlab.ui.control.Label
        UIAxes_2                       matlab.ui.control.UIAxes
        UIAxes                         matlab.ui.control.UIAxes
        AnalyseTab                     matlab.ui.container.Tab
        EditField_2                    matlab.ui.control.EditField
        UITable                        matlab.ui.control.Table
        FileinspectionPanel            matlab.ui.container.Panel
        GetfilesButton                 matlab.ui.control.Button
        PhimaxEditField                matlab.ui.control.NumericEditField
        PhimaxEditFieldLabel           matlab.ui.control.Label
        PhiminEditField                matlab.ui.control.NumericEditField
        PhiminEditFieldLabel           matlab.ui.control.Label
        ThetamaxEditField              matlab.ui.control.NumericEditField
        ThetamaxEditFieldLabel         matlab.ui.control.Label
        ThetaminEditField              matlab.ui.control.NumericEditField
        ThetaminEditFieldLabel         matlab.ui.control.Label
        AngleanalysisPanel             matlab.ui.container.Panel
        LoadButton_2                   matlab.ui.control.Button
        EditField_3                    matlab.ui.control.NumericEditField
        ReturnanglesButton             matlab.ui.control.Button
        Phiofpoint2EditField           matlab.ui.control.NumericEditField
        Phiofpoint2EditFieldLabel      matlab.ui.control.Label
        Thetaofpoint2EditField         matlab.ui.control.NumericEditField
        Thetaofpoint2EditFieldLabel    matlab.ui.control.Label
        Phiofpoint1EditField           matlab.ui.control.NumericEditField
        Phiofpoint1EditFieldLabel      matlab.ui.control.Label
        Thetaofpoint1EditField         matlab.ui.control.NumericEditField
        Thetaofpoint1EditFieldLabel    matlab.ui.control.Label
        UIAxes_3                       matlab.ui.control.UIAxes
        AdvancedanalysisTab            matlab.ui.container.Tab
        GrainclassificationPanel       matlab.ui.container.Panel
        EditField_5                    matlab.ui.control.NumericEditField
        BasisnumberEditField           matlab.ui.control.NumericEditField
        BasisnumberEditFieldLabel      matlab.ui.control.Label
        GrainanalysisButton            matlab.ui.control.Button
        TheoreticalanglesEditField     matlab.ui.control.NumericEditField
        TheoreticalanglesEditFieldLabel  matlab.ui.control.Label
        DeleteButton                   matlab.ui.control.Button
        RownumberEditField             matlab.ui.control.NumericEditField
        RownumberEditFieldLabel        matlab.ui.control.Label
        UITable2                       matlab.ui.control.Table
        DetectpeaksButton              matlab.ui.control.Button
        MosaicrangeEditField           matlab.ui.control.NumericEditField
        MosaicrangeEditFieldLabel      matlab.ui.control.Label
        ThresholdEditField             matlab.ui.control.NumericEditField
        ThresholdEditFieldLabel        matlab.ui.control.Label
        LoadandPlotPanel               matlab.ui.container.Panel
        PhiSlider                      matlab.ui.control.Slider
        PhiLabel                       matlab.ui.control.Label
        PhiEditField                   matlab.ui.control.NumericEditField
        ThetaEditField                 matlab.ui.control.NumericEditField
        ThetaSlider                    matlab.ui.control.Slider
        ThetaSliderLabel               matlab.ui.control.Label
        IntensitycutonDropDown         matlab.ui.control.DropDown
        IntensitycutonDropDownLabel    matlab.ui.control.Label
        EditField_4                    matlab.ui.control.EditField
        LoadButton                     matlab.ui.control.Button
        UIAxes3                        matlab.ui.control.UIAxes
        UIAxes2                        matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        address; %address of current folder
        n_delete;
        starting;
        address_analysis; % Description
        

        adv_n; 
        adv_n2;
        adv_h1;
        adv_h2;
        adv_h3;
        adv_h4;

        
        %To see how many times one has inspect the percentage of grains
        n_grains;
        adv_handles;
        
        adv_grain_information;
        


        ending;

        
        dmax;
        dmin;
        tube;
        pixel;
        n2s;
        resolution;
        grain;
        offset;
        
        adv_step;
        adv_data;      
        
        %Axes range for plotting
        %x_lim_max;
        %x_lim_min;
        %y_lim_max;
        %y_lim_min;
    end
    
    methods (Access = private)
        
        %Convert the cartesian coordinates on the unit sphere to the
        %corresponding (theta, phi) of the spherical polar coordinates, 
        %unit in radians
        
        function poleangles = cartesian2sphere(app,xyz)
            
            theta = asin( sqrt(xyz(1)^2 + xyz(2)^2) / norm(xyz) );
    
            if xyz(3) <= 0
                
                theta = pi - theta;
            end
            
            phi = acos(xyz(1) / sqrt(xyz(1)^2 + xyz(2)^2));
            
            if xyz(2)<0
                
                phi = 2*pi - phi;
                
            end
            
            poleangles = [theta phi];
            
        end
        
        %convert a set of spherical polar coordinates to cartesian
        %coordinates. Phi = polar, Theta = azimuthal, in radian
        function xyz = sphere2cartesian(app,polar,azimuthal)
            
            x = sin(azimuthal) * cos(polar);
    
            y = sin(azimuthal) * sin(polar);
            
            z = cos(azimuthal);
            
            xyz = [x y z]';
            
        end
        
        %Given two phi and theta, find its equivalent cartesian coordinates
        %on the unit sphere, in degrees
        
        function xyz_degree = sphere2cartesian_degree(app,polar,azimuthal)

            x = sind(azimuthal) * cosd(polar);
            
            y = sind(azimuthal) * sind(polar);
            
            z = cosd(azimuthal);
            
            xyz_degree = [x y z]';
        end     

        
        
        %This function is designed to plot the mosaic angle when setting to
        %plot intensity cut on 'Theta'
        
        %To understand how this function works, please refer to
        %grainmaker.m to see how the grain is generated. 
        
        function [delta_phi,counts_cut] = deltaphitheta(app,inputtheta,inputphi)
            
            %Create two lists to store the delta_phi and its
            %corresponding counts
            
            delta_phi = [];
            counts_cut = [];
            
            
            input_theta = inputtheta;
            input_phi = inputphi;
            
            column_number = 0;
            
            %Load step size of grain and each point's neutron counts
            step = app.adv_step;
            data = app.adv_data;
            
            
            %First, consider points with a smaller phi
            for ii = (input_phi/step) : -1 : 0
                column_number = column_number + 1;
                
                delta_phi = [delta_phi -ii*step];
                
                %If phi=0, then it is in the first column and we consider
                %it differently because the number of points in the first
                %column is differet
                if ii == input_phi/step
                    counts_cut = [counts_cut data(input_theta/step + 1)];
                    
                else
                    %Else, find the corresponding counts and append it in
                    %'counts_cut'
                    counts_index = (180/step + 1) + (column_number - 2) * (180/step - 1) + input_theta/step;
                    counts_cut = [counts_cut data(counts_index)];
                end
            end
            
            %The number of remaining points to consider
            number_points = (360/step - 1) - (input_phi/step + 1);
            
            %Loop over the rest points and calculate the delta_phi and
            %neutron counts. 
            for kk = 1 : number_points
                column_number = column_number + 1;
                
                delta_phi = [delta_phi step*kk];
                
                counts_index = (180/step + 1) + (column_number - 2) * (180/step - 1) + input_theta/step;
                counts_cut = [counts_cut data(counts_index)];
            end
            
        end
        
        
        %This function is designed to plot the mosaic angle when setting to
        %plot intensity cut on 'Phi'
        
        %To understand how this function works, please refer to
        %grainmaker.m to see how the grain is generated. 
        
        function [delta_theta,counts_cut] = deltaphitheta2(app,inputtheta,inputphi)
                
                %Create two lists to store the delta_theta and its
                %corresponding counts
                delta_theta = [];
                counts_cut = [];
                
                input_theta = inputtheta;
                input_phi = inputphi;
                
                %Load step size of grain and each point's neutron counts
                step = app.adv_step;
                data = app.adv_data;
                
                
                %If the input phi equals 0, we consider it seperately.
                %Because when we generate the grain, the number of points
                %with phi = 0 is different with the number of points of a
                %different phi value
                if input_phi == 0
                    
                    counts_index = 0;
                    %First, consider all the points with a smaller theta
                    %value than the intersection point of the two straight
                    %lines
                    for ii = (input_theta/step) : -1 : 0
                        
                        counts_index = counts_index + 1;
                        
                        delta_theta = [delta_theta -ii*step];
                        
                        counts_cut = [counts_cut data(counts_index)];
                    end
                    
                    %Calcula the number of the points with phi=0 that are
                    %on the right of the intersection point
                    % (number of all points — the points covered in the previous loop）
                    n_points = (180/step + 1) - (input_theta/step + 1);
                    
                    %Then, loop over the rest points to find the
                    %delta_theta and corresponding neutron counts
                    for jj = 1 : n_points
                       
                        counts_index = counts_index + 1;
                        
                        delta_theta = [delta_theta jj*step];
                        
                        counts_cut = [counts_cut data(counts_index)];
                    end
                    
                else
                    
                    %calculate which column of grain this point corresponds
                    %to
                    column = input_phi/step + 1; 
                    
                    %First, consider the points with a smaller theta
                    counts_index = 0;
                    for ii = (input_theta/step) : -1 : 0
                        
                        counts_index = counts_index + 1;
                        delta_theta = [delta_theta -ii*step];
                        
                        %If theta equals zero, then this point is actually
                        %in the first column
                        if ii == input_theta/step
                            counts_cut = [counts_cut data(counts_index)];
                            
                        else
                            %if not, find its index
                            real_index = (180/step + 1) + (column - 2) * (180/step - 1) + counts_index;
                            %and append the neutron counts to counts_cut
                            counts_cut = [counts_cut data(real_index)];
                        end
                
                    end
                    
                    
                    n_remains = (180/step - 1) - input_theta/step;
                    %from the intersection point, loop to the last point of
                    %the column it belonds
                    for jj = 1 : n_remains
                        
                        delta_theta = [delta_theta jj*step];
                        
                        real_index = real_index + 1;
                        
                        counts_cut = [counts_cut data(real_index)];
                    end
                    
                    %Because the last point of this column is in the first
                    %column
                    delta_theta = [delta_theta (jj+1)*step];
                    counts_cut = [counts_cut data(180/step+1)];    
                end
        end
        
        
        %Draw a circle with radius r centred on (x,y)
        function h = circle(app,x,y,r,ax)
            th = 0:pi/50:2*pi;
            xunit = r * cos(th) + x;
            yunit = r * sin(th) + y;
            h = plot(ax,xunit, yunit)
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: PlotButton
        function PlotButtonPushed2(app, event)
            
            %Distance from neutron source to crystal
            app.n2s= app.DistancebetweenprotonandsamplemEditField.Value;
            
            %the chosen range (dmin, dmax) of the interplanar distance
            app.dmin = app.PlotdfromEditField.Value * 1e-10;
            app.dmax = app.ToEditField_2.Value * 1e-10;
            
            %Number of tubes and pixels per tube
            app.tube = app.NumberoftubesEditField.Value;
            app.pixel = app.PixelspertubeEditField.Value;
            
            %Resolution of visulization
            app.resolution = app.ResolutionDropDown.Value;
            
            %range of the index of pixels to plot (starting, ending)
            app.starting = app.PlotfilesfromEditField.Value;
            app.ending = app.ToEditField.Value;
            
            %Offset angle on Rrot
            app.offset = app.OffsetangleEditField.Value;
            
            %number of pixels to neglect at the ends of each tube,set it to
            %be an integer
            app.n_delete = round(app.NumberofpixelstoneglectattubesendEditField.Value);
            

            
            
            %Set the grain
            grain_resolution = num2str(app.resolution); 
            grain_filename = ['grain(',grain_resolution,').mat'];
            app.grain = load(grain_filename).grain;
            
            %Create nearest neighbour searcher object
            Mdl2 = createns(app.grain);
            
            
            %Mass of neutron and Plank's constant, in standard unit
            M_n = 1.674927471*10^(-27);
            h = 6.62607004*10^(-34);
            
            
            
            
            %Those are the same for all the file
            first_filename = ['ALF',int2str(app.PlotfilesfromEditField.Value),'.nxs'];
            
            distance = h5read(first_filename,'/raw_data_1/instrument/detector_1/distance');
            polarangle = h5read(first_filename,'/raw_data_1/instrument/detector_1/polar_angle');
            aziangle = h5read(first_filename,'/raw_data_1/instrument/detector_1/azimuthal_angle');
            flighttime_pre = h5read(first_filename,'/raw_data_1/instrument/detector_1/time_of_flight');
            
            
            % Convert flighttime from bins to average
            for n=1:length(flighttime_pre)-1
                flighttime(n)=(1/2)*(flighttime_pre(n)+flighttime_pre(n+1));
            end
            flighttime=flighttime';
            
            

            
            %Distance from neutron to detector
            n2d = distance + app.n2s;
            
            %To get the scattering angle theta in radian
            polarangle = polarangle / 2;
            
            
            %Preprocessing the data, neglect any discontinuous point in
            %'aziangle'
            for ii = 2:length(aziangle)-1
                
                if (abs(abs(aziangle(ii)) - abs(aziangle(ii+1))) >= 30) && (abs(abs(aziangle(ii)) - abs(aziangle(ii-1))) >= 30)
                    
                    aziangle(ii) = NaN;
                    
                end
            end
            
            
            %Perform a sanity check
            if app.dmin > app.dmax
                msgbox({'Error: d max needs to be larger than d min'});
            
            elseif app.n_delete < 0
                msgbox({'Error: Number of pixels to neglect can not be negative'});
                
            else

                %Calculate and store the time range accordingto to the d space
                %range
                tmax = ones(app.pixel*app.tube,1);
                tmin = ones(app.pixel*app.tube,1);
                
                for ii = 1:length(aziangle)
                
                    tmax(ii) = (2*sind(polarangle(ii)) * n2d(ii) * M_n * app.dmax) / h;
                    tmin(ii) = (2*sind(polarangle(ii)) * n2d(ii) * M_n * app.dmin) / h;
                    
                end
                
                %Set units in microseconds to compare with dataset
                tmin = tmin * 1000000;
                tmax = tmax * 1000000;
                
                
                %totalcounts is a matrix to store the neutron counts for
                %each corresponding grain's position
                totalcounts = zeros(length(app.grain),1);
                
                %number_adds is a matrix to record, for each grain's spot,
                %how many times the GUI has added a value to it. This is
                %set for average later to prevent overcounting. 
                number_adds = zeros(length(app.grain),1);
                
                %Calculate the range of files to plot
                range = app.ending - app.starting;
                
                %For each file, keep a record of its range in phi and theta
                %those information is stored in matrix 'record'
                
                record = zeros(range,5);
                %First column: index of documents.
                %Second column: maximum of phi.
                %Third column: minimum of phi
                %Forth column: maximum of theta
                %Fifth column: minimum of theta
                
                
                %This is the vector to map as proved in our paper
                initial = [0;1;0];
                
                %Set a progress bar
                progression = waitbar(0,'Ploting in progress');
                
                %Loop over each file
                for kk = app.starting:app.ending
                    
                    waitbar( (kk+1-app.starting)/range,progression)
                    
                    %First, keep a record of the file index
                    record(kk+1-app.starting,1) = kk;
                    
                    %In cases where we have corrupted file, we can neglect
                    %them here
                    if kk == 82562
                        for qq = 2:5
                            record(kk+1-app.starting,qq) = NaN;
                        end
                        continue
                    end
                    
                    
                    %Otherwise, load the experiments data
                    string = int2str(kk);
                    input = ['ALF',string,'.nxs'];
                
                    counts = h5read(input,'/raw_data_1/instrument/detector_1/counts');
                    proton_charge = h5read(input,'/raw_data_1/proton_charge');
                    
                    
                    %Find the intensity detected by each pixel
                    intensity = ones(app.pixel*app.tube,1);
                    for mm = 1:length(aziangle)
                
                        %Find the time channels within the range [tmin, tmax]
                        index = find(flighttime > tmin(mm) & flighttime < tmax(mm));
                
                        target = counts(:,mm);
                        intensity(mm) = sum(target(index));
                    end
                    %Average them by the proton charge
                    intensity = intensity / proton_charge;
                
                    
                    
                    %Load the rotational angle
                    Rlower = h5read(input,'/raw_data_1/selog/Rlower/value_log/value');
                    Rupper = h5read(input,'/raw_data_1/selog/Rupper/value_log/value');
                    Rrot = -h5read(input,'/raw_data_1/selog/Rrot/value_log/value');
                
                    
                    %Take the last value
                    Rrot = Rrot(length(Rrot)) + app.offset;   
                    Rlower = Rlower(length(Rlower));    
                    Rupper = Rupper(length(Rupper));
                    
                    
                
                    %Create rotational matrices
                    inv_Rrot = [cosd(Rrot) -sind(Rrot) 0;
                                    sind(Rrot) cosd(Rrot) 0;
                                    0 0 1];
                    
                    inv_Rlower = [cosd(Rlower) 0 -sind(Rlower);
                                                0 1 0;
                                    sind(Rlower) 0 cosd(Rlower)];
                              
                   
                    inv_Rupper = [1 0 0;
                                     0 cosd(Rupper) -sind(Rupper);
                                     0 sind(Rupper) cosd(Rupper)];
                    
                    
                                 
                    %This matrix is set to record the phi and theta coordinates of each
                    %pixel in the lab frame
                    
                    angle_range = zeros(length(intensity),2);
                    %First column: phi
                    %Second column: theta
                    
                    %Loop over every pixel
                    for ii = 1 : length(intensity)
                        
                        %If this pixel should be neglected, neglect it and
                        %set its angle rangle be NaN
                        if (mod(ii,app.pixel) <= app.NumberofpixelstoneglectattubesendEditField.Value)...
                                || (mod(ii,app.pixel) >= (app.pixel - app.NumberofpixelstoneglectattubesendEditField.Value))
                            
                            for ll = 1:2
                                    
                                angle_range(ii,ll) = NaN;
                            end
                            
                        else
                            
                            %Otherwise, let the five rotational matrices
                            %act on the initial vector to map every pixel
                            %to the lab frame
                            inv_polar = [cosd(polarangle(ii)) sind(polarangle(ii)) 0;
                                          -sind(polarangle(ii)) cosd(polarangle(ii)) 0;
                                          0 0 1];
                

                            inv_azi = [1 0 0;
                                        0 cosd(aziangle(ii)) sind(aziangle(ii));
                                        0 -sind(aziangle(ii)) cosd(aziangle(ii))];
                                    
                
                            pixel_cartesian = inv_Rupper * inv_Rlower * inv_azi * inv_Rrot * inv_polar * initial;
                
                            
                            %Because we have preprocessed the data to set
                            %some of the aziangles to be NaN, if we
                            %encounter such points we neglect them
                            if isnan(pixel_cartesian)
                                
                                
                                for ll = 1:2
                                    
                                    angle_range(ii,ll) = NaN;
                                end
                                
                            else
                                %Find the index of the grain that is
                                %closest to the pixel's cartesian
                                %coordinates
                                Idx = knnsearch(Mdl2, pixel_cartesian');
                                
                                %Add the intensity to the totalcounts and
                                %add one to the number of adds
                                totalcounts(Idx) = totalcounts(Idx) + intensity(ii);
                                number_adds(Idx) = number_adds(Idx) + 1;
                                
                                %convert the mapped pixel's cartesian
                                %coordinates into spherical polar
                                %coordinates, theta and phi. Then, assign
                                %them to the corresponding location in
                                %matrix 'angle_range'
                                corresponding_angle = cartesian2sphere(app,app.grain(Idx,:)) * 180 / pi; 
                                angle_range(ii,1) = corresponding_angle(2);
                                angle_range(ii,2) = corresponding_angle(1);
                                
                            end
                
                        end
                
                    end
                    
                    %Find the range of theta and phi in a given file and
                    %assign the value in record.
                    record(kk+1-app.starting,2) = max(angle_range(:,1));
                    record(kk+1-app.starting,3) = min(angle_range(:,1));
                    record(kk+1-app.starting,4) = max(angle_range(:,2));
                    record(kk+1-app.starting,5) = min(angle_range(:,2));
                end
                
                %save the information
                save('Range_record.mat','record')
                
                
                
                %Average the totalcounts over the number of adds to avoid
                %overcounting
                for mm = 1:length(number_adds)
        
                    if number_adds(mm) ~= 0
                        
                        number_adds(mm) = 1/number_adds(mm);
                    end
                end
                totalcounts = totalcounts .* number_adds;
                
                
                
                %when one point in the grain records zero neutron
                %intensity, its value is set to NaN so it won't show in our
                %plot
                for jj = 1:length(totalcounts)
                    if totalcounts(jj) == 0
                        totalcounts(jj) = NaN;
                    end
                end
                
                %save neutron counts information
                save('totalcounts.mat','totalcounts')
                
                %convert each grain's points' cartesian coordinates into
                %spherical polar coordinates. Store them in matrice
                %'ppangle' and 'aaangle' for plotting. 
                
                ppangle = zeros(length(number_adds),1);
                aaangle = zeros(length(number_adds),1);
                for nn = 1:length(app.grain)
                    
                    angles = cartesian2sphere(app,app.grain(nn,:)) * 180 / pi;
                    
                    ppangle(nn) = angles(2);
                    
                    aaangle(nn) = angles(1);
                end
                
                %Set the range of axes's value
                
                x_lim_max = max(record(:,2)) + 10;
                x_lim_min = min(record(:,3)) - 10;
                
                y_lim_max = max(record(:,4)) + 10;
                y_lim_min = min(record(:,5)) - 10;
                

                %Plot the figrue
                cla(app.UIAxes)
                scatter(app.UIAxes,ppangle,aaangle,[],totalcounts,'filled')
                colorbar(app.UIAxes)
                xlim(app.UIAxes,[x_lim_min x_lim_max]) 
                ylim(app.UIAxes,[y_lim_min y_lim_max])
                
                %Delete the progression bar
                delete(progression)
                
            end

        end

        % Button pushed function: GeneratePoleFigureButton
        function GeneratePoleFigureButtonPushed(app, event)
            
            %Load the information regarding neutron counts detected in each
            %point on the sphere
            data = load('totalcounts.mat').totalcounts;
            %neutron counts of each point
            counts = data(:,1);
            
            %the x,y,z coordinates of each point
            z = app.grain(:,3);
            y = app.grain(:,2);
            x = app.grain(:,1);
            
            %Create three lists to store the coordinates of the projection on the equatorial plane
            %and its corresponding neutron counts
            x_projected = [];
            y_projected = [];
            count_projected = [];
            
            %Sometimes some points are on the lower semisphere or around x
            %or y axis.Therefore, we need to add a correction rotational
            %operator to mapp the points on the sphere before performing 
            %Stereographic projection
            
            correctionx = app.RotatearoundXdegreesEditField.Value;
            correctiony = app.RotatearoundYdegreesEditField.Value;
            
            
                    
            rotation_x = [1 0 0;
                          0 cosd(correctionx) -sind(correctionx);
                          0 sind(correctionx) cosd(correctionx)];
                      
                      
            rotation_y = [cosd(correctiony) 0 -sind(correctiony);
                            0 1 0;
                            sind(correctiony) 0 cosd(correctiony)];
                      
            %test_x, test_y and test_z stores the x,y and z coordinates of
            %points after being rotated by either of the correction matrix
                      
            l = length(x);
            test_x = zeros(l,1);
            test_y = zeros(l,1);
            test_z = zeros(l,1);
            for ii = 1:l
                
                
                
                %Determine from which pole to project here
                
                xyz = (rotation_x * rotation_y * [x(ii); y(ii); z(ii)])';
                
                test_x(ii) = xyz(1);
                test_y(ii) = xyz(2);
                test_z(ii) = xyz(3);
                
                %Neglect all the points in the lower semisphere to avoid
                %overcounting
                if test_z(ii) < 0
                    continue
                else
                    
                    %Find the x and y coordinates of each point's
                    %projection on the equatorial plane
                    spherical_angles = cartesian2sphere(app,xyz);
                    base = sqrt(test_x(ii)^2 + test_y(ii)^2);

                    x_loc = cos(spherical_angles(2)) * base / (abs(test_z(ii)) + 1);
                    y_loc = sin(spherical_angles(2)) * base / (abs(test_z(ii)) + 1);
                    
                    %Store the counts and x,y coordinates of the projection
                    count_projected = [count_projected counts(ii)];
                    x_projected = [x_projected x_loc];
                    y_projected = [y_projected y_loc];
            
                end
            end
            
            %Plot the projection
            cla(app.UIAxes_2);
            scatter(app.UIAxes_2,x_projected,y_projected,8,count_projected,'filled');
            colorbar(app.UIAxes_2)
            xlim(app.UIAxes_2,[-1.2 1.2])
            ylim(app.UIAxes_2,[-1.2 1.2])
            
            %Plot a circle with unit radius
            hold(app.UIAxes_2, 'on' )
            th = 0:pi/50:2*pi;
            xunit = cos(th);
            yunit = sin(th);
            plot(app.UIAxes_2, xunit, yunit);
            hold(app.UIAxes_2, 'on' )
            %Plot the circle's centre
            scatter(app.UIAxes_2,0,0,40,'r','filled');
           
        end

        % Button pushed function: FolderButton
        function FolderButtonPushed(app, event)
            %Add the required path
            app.address = uigetdir;
            addpath(app.address)
            app.EditField.Value = convertCharsToStrings(app.address);
        end

        % Button pushed function: LoadButton_2
        function LoadButton_2Pushed(app, event)

            %Add the required path
            third_address = uigetdir;
            addpath(third_address)
            app.EditField_2.Value = convertCharsToStrings(third_address);
           
            %Load information about neutron counts
            data_analysis = load('totalcounts.mat').totalcounts;
            
            %Determine which grain we were using when generating the
            %'totalcounts.mat' and load the grain
            if length(data_analysis) == 16022
                grain_analysis = load('grain(2).mat').grain;
            else
                grain_analysis = load('grain(1).mat').grain;
            end
            
            %Convert the cartesian coordinates of each point in the grain
            %into spherical polar coordinates       
            ppangle = zeros(length(grain_analysis),1);
            aaangle = zeros(length(grain_analysis),1);
            for nn = 1:length(grain_analysis)
                
                angles = cartesian2sphere(app,grain_analysis(nn,:)) * 180 / pi;
                
                ppangle(nn) = angles(2);
                
                aaangle(nn) = angles(1);
            end
            
            %load the records of the range of angles
            records = load('Range_record.mat').record;
            
            %The range of axes' value
            x_lim_max = max(records(:,2)) + 10;
            x_lim_min = min(records(:,3)) - 10;
                
            y_lim_max = max(records(:,4)) + 10;
            y_lim_min = min(records(:,5)) - 10;
            

            
            %Plot the neutron counts information like what we did in the
            %second tab
            cla(app.UIAxes_3)
            scatter(app.UIAxes_3,ppangle,aaangle,[],data_analysis,'filled')
            
            xlim(app.UIAxes_3,[x_lim_min x_lim_max]) 
            ylim(app.UIAxes_3,[y_lim_min y_lim_max])
            
            colorbar(app.UIAxes_3)
        end

        % Button pushed function: ReturnanglesButton
        function ReturnanglesButtonPushed(app, event)
            %a,b,c,d are the theta and phi value of two points of interest
            a = app.Thetaofpoint1EditField.Value;
            b = app.Phiofpoint1EditField.Value;
            
            c = app.Thetaofpoint2EditField.Value;
            d = app.Phiofpoint2EditField.Value;
            
            
            %Calculate the cartesian coordinates of the two points
            input1 = [a b];
            input2 = [c d];
            
            input1_car = sphere2cartesian_degree(app,input1(1),input1(2));
            input2_car = sphere2cartesian_degree(app,input2(1),input2(2));
            
            %Calculate the angle between them and print it on the GUI
            app.EditField_3.Value = acosd(dot(input1_car,input2_car)/norm(input1_car)/norm(input2_car));
        end

        % Button pushed function: GetfilesButton
        function GetfilesButtonPushed(app, event)
            %load the records of the range of angles
            records = load('Range_record.mat').record;
            
            
            %Read the maximum and minimum value of theta, phi
            theta_min = app.ThetaminEditField.Value;
            theta_max = app.ThetamaxEditField.Value;
            
            phi_min = app.PhiminEditField.Value;
            phi_max = app.PhimaxEditField.Value;
            
            
            %Perform sanity checks
            if theta_min > theta_max
                msgbox({'Error: theta_max needs to be larger than theta_min'});
            
            elseif phi_min > phi_max
                msgbox({'Error: phi_max needs to be larger than phi_min'});
                
            else
                
                %Only select files that contains points within this range,
                %and append the file index to 'filtered_list'
                
                filtered_list = [];
                for ii = 1:length(records)
                    
                    if phi_max >= records(ii,3) && phi_min <= records(ii,2)...
                            && theta_max >= records(ii,5) && theta_min <= records(ii,4)
                        
                        filtered_list = [filtered_list;records(ii,:)];
                    end
                end
                
                %Clear the table and print the list
                app.UITable.Data(:,:,:,:,:); 
                
                if isempty(filtered_list)
                    msgbox({'No files found within this range'});
                else
                    app.UITable.Data = filtered_list;
                end
            end
        end

        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
            %Add the correct path
            fourth_address = uigetdir;
            addpath(fourth_address)
            app.EditField_4.Value = convertCharsToStrings(fourth_address);
           
            %Determine which grain we were using when generating
            %'totalcounts.mat'
            app.adv_data = load('totalcounts.mat').totalcounts;
            if length(app.adv_data) == 16022
                adv_grain = load('grain(2).mat').grain;
                app.adv_step = 2;
            else
                adv_grain = load('grain(1).mat').grain;
                app.adv_step = 1;
            end
            
            %Convert the grain's cartesian coordinates to spherical polar
            %coordinates
            adv_ppangle = zeros(length(adv_grain),1);
            adv_aaangle = zeros(length(adv_grain),1);

            for nn = 1:length(adv_grain)
                
                angles = cartesian2sphere(app,adv_grain(nn,:)) * 180 / pi;
                
                adv_ppangle(nn) = angles(2);
                
                adv_aaangle(nn) = angles(1);
            end
            
            
            %load the records of the range of angles
            records = load('Range_record.mat').record;
            
            %The range of axes' value
            x_lim_max = max(records(:,2)) + 10;
            x_lim_min = min(records(:,3)) - 10;
                
            y_lim_max = max(records(:,4)) + 10;
            y_lim_min = min(records(:,5)) - 10;
            

            %Plot the neutron intensity on the sphere
            cla(app.UIAxes2)
            scatter(app.UIAxes2,adv_ppangle,adv_aaangle,[],app.adv_data,'filled')
            xlim(app.UIAxes2,[x_lim_min x_lim_max])
            ylim(app.UIAxes2,[y_lim_min y_lim_max])
            colorbar(app.UIAxes2)
            hold(app.UIAxes2,'on')
            
            %To record the number of times we slide Phi
            app.adv_n = 0;
            %To record the number of times we slide Theta
            app.adv_n2 = 0;
            
            app.n_grains = 0;
            app.adv_handles = [];
            
            
        end

        % Value changing function: PhiSlider
        function PhiSliderValueChanging(app, event)
            
            changingValue = round(event.Value);
            
            %Only select the multiple of 2 as the phi value, if two
            %adjacent points on the grain are 2 degrees apart
            if mod(changingValue,2) == 1 && app.adv_step == 2
                changingValue = changingValue - 1; 
            end
            
            %print the phi value
            app.PhiEditField.Value = changingValue;
            
            %update the number of times we have slided phi
            app.adv_n = app.adv_n + 1;
            
            %Two straight lines for plot, 
            %theta = app.ThetaEditField.Value and phi = app.PhiEditField.Value
            theta = ones(1000,1) * app.ThetaEditField.Value;
            phi = ones(1000,1) * app.PhiEditField.Value;
            
            %Set the domain range of theta and phi
            theta_x = linspace(0,180,1000);
            phi_x = linspace(0,360,1000);
            
            
            
            if app.adv_n == 1
                %If this is the first time we slide phi but not the first
                %time we slide Theta, delete the two straight lines from
                %previous plot
                if app.adv_n2 ~= 1
                    delete(app.adv_h3)
                    delete(app.adv_h4)
                end
                %and then plot the two straight lines
                app.adv_h1 = plot(app.UIAxes2,phi_x,theta,'r');
                hold(app.UIAxes2,'on')
                app.adv_h2 = plot(app.UIAxes2,phi,theta_x,'r');
            else
                
                %else if it is not the first time we move both slider,
                %delete all the previous lines
                if app.adv_n2 ~= 1
                    delete(app.adv_h3)
                    delete(app.adv_h4)
                end
                delete(app.adv_h1)
                delete(app.adv_h2)
                
                %Then plot the new two lines
                app.adv_h1 = plot(app.UIAxes2,phi_x,theta,'r');
                hold(app.UIAxes2,'on')
                app.adv_h2 = plot(app.UIAxes2,phi,theta_x,'r');
            end
            
            %Generate the x,y data for plotting depending on whether we
            %choose to plot intensity cut on 'Theta' or 'Phi'
            if strcmp(app.IntensitycutonDropDown.Value,'Theta') == 1

                [x,y] = deltaphitheta(app,app.ThetaEditField.Value,app.PhiEditField.Value);
                %xlim(app.UIAxes3,[-360,360])
                xlabel(app.UIAxes3,'\Delta{\Phi}')
                
            else
                [x,y] = deltaphitheta2(app,app.ThetaEditField.Value,app.PhiEditField.Value);
                %xlim(app.UIAxes3,[-180,180])
                xlabel(app.UIAxes3,'\Delta{\Theta}')
            end
    
            %Clear the figure and plot
            cla(app.UIAxes3)
            plot(app.UIAxes3,x,y)
        end

        % Value changing function: ThetaSlider
        function ThetaSliderValueChanging(app, event)
            changingValue = round(event.Value);
            
            %Only select the multiple of 2 as the phi value, if two
            %adjacent points on the grain are 2 degrees apart
            if mod(changingValue,2) == 1 && app.adv_step == 2
                changingValue = changingValue - 1; 
            end

            app.ThetaEditField.Value = changingValue;
            
            %update the number of times we have slided theta
            app.adv_n2 = app.adv_n2 + 1;
            
            %Two straight lines for plot, 
            %theta = app.ThetaEditField.Value and phi = app.PhiEditField.Value
            theta = ones(1000,1) * app.ThetaEditField.Value;
            phi = ones(1000,1) * app.PhiEditField.Value;
            
            %Set the domain range of theta and phi
            theta_x = linspace(0,180,1000);
            phi_x = linspace(0,360,1000);
            
            
            if app.adv_n2 == 1
                %If this is the first time we slide Theta but not the first
                %time we slide Phi
                if app.adv_n ~= 1
                    %Delete the two lines from previous plot
                    delete(app.adv_h1)
                    delete(app.adv_h2)
                end
                
                %Plot the new two lines
                app.adv_h3 = plot(app.UIAxes2,phi_x,theta,'r');
                hold(app.UIAxes2,'on')
                app.adv_h4 = plot(app.UIAxes2,phi,theta_x,'r');
            else
                
                %If this is not the first time we slide theta and phi,
                %delete all previous plotted straight lines
                if app.adv_n ~= 1
                    delete(app.adv_h1)
                    delete(app.adv_h2)
                end
                delete(app.adv_h3)
                delete(app.adv_h4)
                
                %plot two new lines
                app.adv_h3 = plot(app.UIAxes2,phi_x,theta,'r');
                hold(app.UIAxes2,'on')
                app.adv_h4 = plot(app.UIAxes2,phi,theta_x,'r');
            end
            
            %Generate the x,y data for plotting depending on whether we
            %choose to plot intensity cut on 'Theta' or 'Phi'
            if strcmp(app.IntensitycutonDropDown.Value,'Theta') == 1

                [x,y] = deltaphitheta(app,app.ThetaEditField.Value,app.PhiEditField.Value);
                %xlim(app.UIAxes3,[-360,360])
                xlabel(app.UIAxes3,'\Delta{\Phi}')
                
            else
                [x,y] = deltaphitheta2(app,app.ThetaEditField.Value,app.PhiEditField.Value);
                %xlim(app.UIAxes3,[-180,180])
                xlabel(app.UIAxes3,'\Delta{\Theta}')
                
            end
            
            %Clear the figure and plot
            cla(app.UIAxes3)
            plot(app.UIAxes3,x,y)
            
        end

        % Button pushed function: PlotButton_2
        function PlotButton_2Pushed(app, event)
            
            %The range of interplanar distance is converted to a unit of
            %10^(-10) m
            sin_dmax = app.dmaxEditField.Value * 10^(-10);
            sin_dmin = app.dminEditField.Value * 10^(-10);
            
            
            %Perform a sanity check
            if sin_dmin > sin_dmax
                msgbox({'Error: d max needs to be larger than d min'});
            else
        
                %Mass of neutron and Plank's constant, in standard unit
                M_n = 1.674927471*10^(-27);
                h = 6.62607004*10^(-34);
                
                %Initiate some instruments' parameters
                sin_tube = app.NumberoftubesEditField_2.Value;
                sin_pixel = app.PixelspertubeEditField_2.Value;
                
                %Distance from neutron to sample. unit in meters
                sin_n2s = app.DistancebetweenprotonandsamplemEditField_2.Value;
                
                
                %Read the information from the input file's index
                file_index = app.FileindexEditField.Value;
                single_input = ['ALF',num2str(file_index),'.nxs'];
                
                distance = h5read(single_input,'/raw_data_1/instrument/detector_1/distance');
                polarangle = h5read(single_input,'/raw_data_1/instrument/detector_1/polar_angle');
                azimuthalangle = h5read(single_input,'/raw_data_1/instrument/detector_1/azimuthal_angle');
                flighttime_pre = h5read(single_input,'/raw_data_1/instrument/detector_1/time_of_flight');
                
                %Neutron counts
                counts = h5read(single_input,'/raw_data_1/instrument/detector_1/counts');
                %Magnitude of proton charge which created the incoming neutrons
                proton_charge = h5read(single_input,'/raw_data_1/proton_charge');
    
                
                % Convert flighttime from bins to average
                for n=1:length(flighttime_pre)-1
                    flighttime(n)=(1/2)*(flighttime_pre(n)+flighttime_pre(n+1));
                end
                flighttime=flighttime';
                
                
                polarangle = polarangle / 2;
                
                %Initiate three matrix with size (n_pixel, n_tube) to store the
                %polar angle, azimuthal angle and distance to the crystal for
                %each pixel
                d_matrix = ones(sin_pixel,sin_tube);
                polar_matrix = ones(sin_pixel,sin_tube);
                azi_matrix = ones(sin_pixel,sin_tube);
                
                %Calculate the distance, polar and azimuthal angle corresponding to each
                %pixel
                
                %Loop over each tube
                for ii = 1:sin_tube
                    %And update each pixel
                    n = 0;
                    for jj = (ii-1) * sin_pixel + 1 : ii * sin_pixel
                        n = n+1;
                        d_matrix(n,ii) = distance(jj);
                        polar_matrix(n,ii) = polarangle(jj);
                        azi_matrix(n,ii) = azimuthalangle(jj);
                    end
                end
                
                %Update 'd_matrix' to represent the distance from source to detector
                d_matrix = d_matrix + sin_n2s;
                
                
                
                %Calculate the range of flying time of detected neutrons of
                %each pixel and store them in matrices tmax and tmin
                tmax = ones(sin_pixel,sin_tube);
                tmin = ones(sin_pixel,sin_tube);
                for ii = 1:sin_tube
                    
                    for jj = 1:sin_pixel
                
                        tmax(jj,ii) = (2*sind(polar_matrix(jj,ii)) * d_matrix(jj,ii) * M_n * sin_dmax) / h;
                        tmin(jj,ii) = (2*sind(polar_matrix(jj,ii)) * d_matrix(jj,ii) * M_n * sin_dmin) / h;
                    end
                end
                
                
                %Set units in microseconds to compare with dataset
                tmin = tmin * 1000000;
                tmax = tmax * 1000000;
                
    
    
                intensity = ones(sin_pixel,sin_tube);   
                
                %Find the time channels within the range [tmin, tmax], and
                %store the sum of neutron counts that arrived within this time range
                %in matrix 'intensity'
                for mm = 1:sin_tube
                    
                    for nn = 1:sin_pixel
                        
                        location = (mm-1) * sin_pixel + nn;
                        
                        target = counts(:,location);
            
                        index = find(flighttime > tmin(nn,mm) & flighttime < tmax(nn,mm));
             
                        intensity(nn,mm) = sum(target(index));
                        
                    end
                end
                
                %Clear the figure, to avoid overlap with previous plots
                cla(app.UIAxes4)
                %Normalize the neutron intensity with the proton charge
                intensity = intensity / proton_charge;
                
                %Plot the detector
                imagesc(app.UIAxes4,intensity)
                xlim(app.UIAxes4,[0,sin_tube])
                ylim(app.UIAxes4,[0,sin_pixel])
                colorbar(app.UIAxes4)
                
                
                %Read the angles information and show them in the GUI
                Rlower = h5read(single_input,'/raw_data_1/selog/Rlower/value_log/value');
                Rupper = h5read(single_input,'/raw_data_1/selog/Rupper/value_log/value');        
                Rrot = h5read(single_input,'/raw_data_1/selog/Rrot/value_log/value');
            
                %Take the last value
                Rrot = double(Rrot(length(Rrot)));
                Rupper = double(Rupper(length(Rupper)));
                Rlower = double(Rlower(length(Rlower)));
    
                %Print them out
                app.RrotEditField.Value = Rrot;   
                app.RlowerEditField.Value = Rlower;    
                app.RupperEditField.Value = Rupper;
            end

        end

        % Button pushed function: DetectpeaksButton
        function DetectpeaksButtonPushed(app, event)
            
            %To keep a record to see how many times we have detected peaks
            app.n_grains = app.n_grains + 1;
            
            %Detect all the peaks above this threshold
            threshold = app.ThresholdEditField.Value;
            
            %width of each peak
            radius = app.MosaicrangeEditField.Value * pi / 180;
            

            grain_data = app.adv_data;
            %Change all the nan values to zero
            for iii = 1 : length(grain_data)
                if isnan(grain_data(iii))
                    grain_data(iii) = 0;
                end
            end
            
            %Determine which grain we were using when visualizing the data
            if length(grain_data) == 16022
                adv_grain = load('grain(2).mat').grain;
            else
                adv_grain = load('grain(1).mat').grain;
            end
            
            maximum = max(grain_data);
            
            %perform a sanity check
            if maximum <= threshold
                msgbox('No peaks found, please try a lower threshold');
            else
            
                %Store the information about grain poitns
                %rownumber; theta; phi; counts
                app.adv_grain_information = [];
                
                rownumber = 0;
                for ii = ceil(maximum) : -1 : threshold 
                    
                    [local_max,location] = max(grain_data);
    
                    if ii <= local_max
                        
                        rownumber = rownumber+1;
                        
                        %Theta phi
                        angles = cartesian2sphere(app,adv_grain(location,:)) * 180 / pi;
                        
                        
                        %find the surrounding points
                        surrounds_location = [];
                        
                        %The counts of peak
                        peak_counts = 0;
                        for jj = 1:length(adv_grain)
                            
                            %Find all information about the peak
                            if norm(adv_grain(jj,:) - adv_grain(location,:)) <= radius
                                surrounds_location = [surrounds_location jj];
                                peak_counts = peak_counts + grain_data(jj);
                            end
                        end
                          
                        app.adv_grain_information = [app.adv_grain_information; rownumber angles(1) angles(2) peak_counts];
                        
                        %Destory the peak by setting its surrounding within the
                        %mosaic range to zero
                        grain_data(surrounds_location) = 0;
                    end
    
                end
                
                %Clear the table
                app.UITable2.Data(:,:,:,:);
                %Print the information of detected peak
                app.UITable2.Data = app.adv_grain_information;
                
                
                %Draw a circle around each detected peak
                th = 0:pi/50:2*pi;
                for kk = 1:length(app.adv_grain_information(:,1))
                    
                    %if this is the first time we detect peaks
                    if app.n_grains == 1
                        
                        hold(app.UIAxes2,'on')
                        
                        xunit = radius * 180 / pi * cos(th) + app.adv_grain_information(kk,2);
                        yunit = radius * 180 / pi * sin(th) + app.adv_grain_information(kk,3);
                        
                        h = plot(app.UIAxes2,yunit, xunit,'r','linewidth',1);
    
                        app.adv_handles = [app.adv_handles h];
                    else
                        %If not, we first delete the previous circles
                        if kk == 1
                            delete(app.adv_handles)
                        end
                        
                        %Then plot new ones
                        hold(app.UIAxes2,'on')
                        xunit = radius * 180 / pi * cos(th) + app.adv_grain_information(kk,2);
                        yunit = radius * 180 / pi * sin(th) + app.adv_grain_information(kk,3);
                        h = plot(app.UIAxes2,yunit, xunit,'r','linewidth',1);
    
                        app.adv_handles(kk) = h;
                    end
                end
            end
        end

        % Button pushed function: DeleteButton
        function DeleteButtonPushed(app, event)
            
            %Choose the row's index to delete
            n_todelete = app.RownumberEditField.Value;
            for ii = 1 : length(app.UITable2.Data(:,1))
                %Delete the line with the corresponding row number
                if app.UITable2.Data(ii,1) == n_todelete
                    app.UITable2.Data(ii,:,:,:) = [];
                    break
                end
                
            end
            %Delete the circles of the deleted peak
            delete(app.adv_handles(n_todelete));
        end

        % Button pushed function: GrainanalysisButton
        function GrainanalysisButtonPushed(app, event)
            
            %Final_info contains the selected peaks
            final_info = app.UITable2.Data;
            %The theoretical angles between peaks of this certain atom
            %planes
            angle_theory = app.TheoreticalanglesEditField.Value;
            
            
            %Find the index of the basis
            basisnumber1 = app.BasisnumberEditField.Value;
            basisnumber2 = app.EditField_5.Value;
            
            %Rownumber, theta, phi, counts
            basis = [];
            %Left overs
            todetermine = [];
            
            %Allocate the information from the table to corresponding lists
            for ii = 1 : length(final_info(:,1))
                if final_info(ii,1) == basisnumber1 || final_info(ii,1) == basisnumber2
                    basis = [basis; final_info(ii,:)];
                else
                    todetermine = [todetermine; final_info(ii,:)];
                end
            end
            
            %Peaks of the same grain needs to be in the same line, first
            %calculate the basis's slope
            standard_slope = atand((basis(2,2) - basis(1,2)) / (basis(2,3) - basis(1,3)));
            
 
            %Store the index of rows that are the same grain
            selected = [];

            
            for ii = 1 : length(todetermine(:,1))
                
                %Only when it satisfies both basis will it be added
                double_check = 0;
                
                %Basis only has two rows
                for jj = 1:2
                    
                    %phi theta of input location
                    input1 = [todetermine(ii,3) todetermine(ii,2)];
                    
                    %phi theta of basis
                    input2 = [basis(jj,3) basis(jj,2)];
                    
                    %Convert them to cartesian coordinates
                    input1_car = sphere2cartesian_degree(app,input1(1),input1(2)); 
                    input2_car = sphere2cartesian_degree(app,input2(1),input2(2));
                    
                    %Calculate the angle difference
                    angle_diff = acosd(dot(input1_car,input2_car)/norm(input1_car)/norm(input2_car));
                    
                    %The slope of input location and the jth peak
                    check_slope = atand((todetermine(ii,2) - basis(jj,2)) / (todetermine(ii,3) - basis(jj,3)));
                    %difference in slope
                    slope_diff = abs(check_slope - standard_slope);
                    
                    
                    %N.B angle_diff is positive by design
                    
                    if mod(angle_diff, angle_theory) <= app.MosaicrangeEditField.Value ...
                            || mod(angle_diff, angle_theory) >= (angle_theory - app.MosaicrangeEditField.Value) ...
                            && slope_diff <= app.MosaicrangeEditField.Value
                        
                        %Make sure it satisfies condition with both bases. 
                        double_check = double_check + 1;
                        
                        if double_check == 2
                            selected = [selected; todetermine(ii,:)];
                        end
                    end
                end
            end
            
            %If no peaks found, return this message
            if isempty(selected)
                msgbox({'No peaks of the same grain are found'});
            else
                
                %Else, return the message containing the index of peaks
                %from the same grain and the average neutron counts of them
                message = [];
                average = 0;
                for mm = 1 : length(selected(:,1))
                    message = [message;num2str(selected(mm,1))];
                    average = average + selected(mm,4);
                end
                
                average = average / length(selected(:,1));
                
                array = ['are also of the same grain with average counts ',num2str(average)];
                
                msgbox({'Peaks with row number:';message; array});
            end
        end

        % Button pushed function: FolderButton_2
        function FolderButton_2Pushed(app, event)
            %Get the address needed
            first_address = uigetdir;
            %Add to path
            addpath(convertCharsToStrings(first_address));
            %Print the address in the GUI
            app.EditField_6.Value = convertCharsToStrings(first_address);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 766 601];
            app.UIFigure.Name = 'UI Figure';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 -19 766 621];

            % Create SinglefilevisualizationTab
            app.SinglefilevisualizationTab = uitab(app.TabGroup);
            app.SinglefilevisualizationTab.Title = 'Single file visualization';

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.SinglefilevisualizationTab);
            title(app.UIAxes4, 'Detector')
            xlabel(app.UIAxes4, 'Tubes')
            ylabel(app.UIAxes4, 'Pixels')
            app.UIAxes4.XTickLabelRotation = 0;
            app.UIAxes4.YTickLabelRotation = 0;
            app.UIAxes4.ZTickLabelRotation = 0;
            app.UIAxes4.Position = [12 193 515 380];

            % Create PlotsettingsPanel
            app.PlotsettingsPanel = uipanel(app.SinglefilevisualizationTab);
            app.PlotsettingsPanel.TitlePosition = 'centertop';
            app.PlotsettingsPanel.Title = 'Plot settings';
            app.PlotsettingsPanel.Position = [538 340 218 198];

            % Create PlotButton_2
            app.PlotButton_2 = uibutton(app.PlotsettingsPanel, 'push');
            app.PlotButton_2.ButtonPushedFcn = createCallbackFcn(app, @PlotButton_2Pushed, true);
            app.PlotButton_2.Position = [70 20 100 22];
            app.PlotButton_2.Text = 'Plot';

            % Create FileindexEditFieldLabel
            app.FileindexEditFieldLabel = uilabel(app.PlotsettingsPanel);
            app.FileindexEditFieldLabel.HorizontalAlignment = 'right';
            app.FileindexEditFieldLabel.Position = [54 142 57 22];
            app.FileindexEditFieldLabel.Text = 'File index';

            % Create FileindexEditField
            app.FileindexEditField = uieditfield(app.PlotsettingsPanel, 'numeric');
            app.FileindexEditField.ValueDisplayFormat = '%.0f';
            app.FileindexEditField.Position = [117 142 44 22];

            % Create dminEditFieldLabel
            app.dminEditFieldLabel = uilabel(app.PlotsettingsPanel);
            app.dminEditFieldLabel.HorizontalAlignment = 'right';
            app.dminEditFieldLabel.Position = [65 104 35 22];
            app.dminEditFieldLabel.Text = 'd min';

            % Create dminEditField
            app.dminEditField = uieditfield(app.PlotsettingsPanel, 'numeric');
            app.dminEditField.Position = [117 104 44 22];

            % Create dmaxEditFieldLabel
            app.dmaxEditFieldLabel = uilabel(app.PlotsettingsPanel);
            app.dmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.dmaxEditFieldLabel.Position = [66 63 38 22];
            app.dmaxEditFieldLabel.Text = 'd max';

            % Create dmaxEditField
            app.dmaxEditField = uieditfield(app.PlotsettingsPanel, 'numeric');
            app.dmaxEditField.Position = [117 63 44 22];

            % Create NumberoftubesEditField_2Label
            app.NumberoftubesEditField_2Label = uilabel(app.SinglefilevisualizationTab);
            app.NumberoftubesEditField_2Label.HorizontalAlignment = 'right';
            app.NumberoftubesEditField_2Label.Position = [245 77 94 22];
            app.NumberoftubesEditField_2Label.Text = 'Number of tubes';

            % Create NumberoftubesEditField_2
            app.NumberoftubesEditField_2 = uieditfield(app.SinglefilevisualizationTab, 'numeric');
            app.NumberoftubesEditField_2.Position = [363 77 58 22];
            app.NumberoftubesEditField_2.Value = 37;

            % Create PixelspertubeEditField_2Label
            app.PixelspertubeEditField_2Label = uilabel(app.SinglefilevisualizationTab);
            app.PixelspertubeEditField_2Label.HorizontalAlignment = 'right';
            app.PixelspertubeEditField_2Label.Position = [254 34 85 22];
            app.PixelspertubeEditField_2Label.Text = 'Pixels per tube';

            % Create PixelspertubeEditField_2
            app.PixelspertubeEditField_2 = uieditfield(app.SinglefilevisualizationTab, 'numeric');
            app.PixelspertubeEditField_2.Position = [363 34 58 22];
            app.PixelspertubeEditField_2.Value = 64;

            % Create DistancebetweenprotonandsamplemEditField_2Label
            app.DistancebetweenprotonandsamplemEditField_2Label = uilabel(app.SinglefilevisualizationTab);
            app.DistancebetweenprotonandsamplemEditField_2Label.Position = [119 121 232 22];
            app.DistancebetweenprotonandsamplemEditField_2Label.Text = 'Distance between proton and sample (m)';

            % Create DistancebetweenprotonandsamplemEditField_2
            app.DistancebetweenprotonandsamplemEditField_2 = uieditfield(app.SinglefilevisualizationTab, 'numeric');
            app.DistancebetweenprotonandsamplemEditField_2.ValueDisplayFormat = '%.2f';
            app.DistancebetweenprotonandsamplemEditField_2.Position = [363 121 58 22];
            app.DistancebetweenprotonandsamplemEditField_2.Value = 14.863;

            % Create RrotEditFieldLabel
            app.RrotEditFieldLabel = uilabel(app.SinglefilevisualizationTab);
            app.RrotEditFieldLabel.HorizontalAlignment = 'right';
            app.RrotEditFieldLabel.Position = [574 270 28 22];
            app.RrotEditFieldLabel.Text = 'Rrot';

            % Create RrotEditField
            app.RrotEditField = uieditfield(app.SinglefilevisualizationTab, 'numeric');
            app.RrotEditField.Position = [617 270 100 22];

            % Create RupperEditFieldLabel
            app.RupperEditFieldLabel = uilabel(app.SinglefilevisualizationTab);
            app.RupperEditFieldLabel.HorizontalAlignment = 'right';
            app.RupperEditFieldLabel.Position = [557 231 45 22];
            app.RupperEditFieldLabel.Text = 'Rupper';

            % Create RupperEditField
            app.RupperEditField = uieditfield(app.SinglefilevisualizationTab, 'numeric');
            app.RupperEditField.Position = [617 231 100 22];

            % Create RlowerEditFieldLabel
            app.RlowerEditFieldLabel = uilabel(app.SinglefilevisualizationTab);
            app.RlowerEditFieldLabel.HorizontalAlignment = 'right';
            app.RlowerEditFieldLabel.Position = [559 193 43 22];
            app.RlowerEditFieldLabel.Text = 'Rlower';

            % Create RlowerEditField
            app.RlowerEditField = uieditfield(app.SinglefilevisualizationTab, 'numeric');
            app.RlowerEditField.Position = [617 193 100 22];

            % Create FolderButton_2
            app.FolderButton_2 = uibutton(app.SinglefilevisualizationTab, 'push');
            app.FolderButton_2.ButtonPushedFcn = createCallbackFcn(app, @FolderButton_2Pushed, true);
            app.FolderButton_2.Position = [119 156 74 22];
            app.FolderButton_2.Text = 'Folder';

            % Create EditField_6
            app.EditField_6 = uieditfield(app.SinglefilevisualizationTab, 'text');
            app.EditField_6.Position = [217 156 204 22];

            % Create TextArea_2
            app.TextArea_2 = uitextarea(app.SinglefilevisualizationTab);
            app.TextArea_2.HorizontalAlignment = 'center';
            app.TextArea_2.Position = [451 34 295 123];
            app.TextArea_2.Value = {'This software is developed by Zihao Liu under the supervision of Prof Christopher Stock at the University of Edinburgh. '; ''; 'To report any bug or provide future suggestions please contact:'; ''; 's1701703@ed.ac.uk'};

            % Create PlotTab
            app.PlotTab = uitab(app.TabGroup);
            app.PlotTab.Title = 'Plot';

            % Create UIAxes
            app.UIAxes = uiaxes(app.PlotTab);
            xlabel(app.UIAxes, '\Phi')
            ylabel(app.UIAxes, {'\Theta'; ''})
            app.UIAxes.XTickLabelRotation = 0;
            app.UIAxes.YTickLabelRotation = 0;
            app.UIAxes.ZTickLabelRotation = 0;
            app.UIAxes.Position = [431 295 301 292];

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.PlotTab);
            title(app.UIAxes_2, 'Pole Figure')
            ylabel(app.UIAxes_2, {''; ''})
            app.UIAxes_2.XTickLabelRotation = 0;
            app.UIAxes_2.YTickLabelRotation = 0;
            app.UIAxes_2.ZTickLabelRotation = 0;
            app.UIAxes_2.Position = [431 34 311 247];

            % Create ExperimentSetupPanel
            app.ExperimentSetupPanel = uipanel(app.PlotTab);
            app.ExperimentSetupPanel.TitlePosition = 'centertop';
            app.ExperimentSetupPanel.Title = 'Experiment Set-up';
            app.ExperimentSetupPanel.Position = [11 428 407 159];

            % Create DistancebetweenprotonandsamplemEditFieldLabel
            app.DistancebetweenprotonandsamplemEditFieldLabel = uilabel(app.ExperimentSetupPanel);
            app.DistancebetweenprotonandsamplemEditFieldLabel.Position = [15 106 232 22];
            app.DistancebetweenprotonandsamplemEditFieldLabel.Text = 'Distance between proton and sample (m)';

            % Create DistancebetweenprotonandsamplemEditField
            app.DistancebetweenprotonandsamplemEditField = uieditfield(app.ExperimentSetupPanel, 'numeric');
            app.DistancebetweenprotonandsamplemEditField.ValueDisplayFormat = '%.2f';
            app.DistancebetweenprotonandsamplemEditField.Position = [256 106 111 22];
            app.DistancebetweenprotonandsamplemEditField.Value = 14.863;

            % Create NumberoftubesEditFieldLabel
            app.NumberoftubesEditFieldLabel = uilabel(app.ExperimentSetupPanel);
            app.NumberoftubesEditFieldLabel.Position = [15 76 232 22];
            app.NumberoftubesEditFieldLabel.Text = 'Number of tubes';

            % Create NumberoftubesEditField
            app.NumberoftubesEditField = uieditfield(app.ExperimentSetupPanel, 'numeric');
            app.NumberoftubesEditField.ValueDisplayFormat = '%.0f';
            app.NumberoftubesEditField.Position = [256 76 111 22];
            app.NumberoftubesEditField.Value = 37;

            % Create PixelspertubeEditFieldLabel
            app.PixelspertubeEditFieldLabel = uilabel(app.ExperimentSetupPanel);
            app.PixelspertubeEditFieldLabel.Position = [14 43 232 22];
            app.PixelspertubeEditFieldLabel.Text = 'Pixels per tube';

            % Create PixelspertubeEditField
            app.PixelspertubeEditField = uieditfield(app.ExperimentSetupPanel, 'numeric');
            app.PixelspertubeEditField.ValueDisplayFormat = '%.0f';
            app.PixelspertubeEditField.Position = [255 43 111 22];
            app.PixelspertubeEditField.Value = 64;

            % Create OffsetangleEditFieldLabel
            app.OffsetangleEditFieldLabel = uilabel(app.ExperimentSetupPanel);
            app.OffsetangleEditFieldLabel.HorizontalAlignment = 'right';
            app.OffsetangleEditFieldLabel.Position = [8 11 72 22];
            app.OffsetangleEditFieldLabel.Text = 'Offset angle';

            % Create OffsetangleEditField
            app.OffsetangleEditField = uieditfield(app.ExperimentSetupPanel, 'numeric');
            app.OffsetangleEditField.Position = [255 11 113 22];

            % Create PlotButton
            app.PlotButton = uibutton(app.PlotTab, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @PlotButtonPushed2, true);
            app.PlotButton.Position = [151 198 100 22];
            app.PlotButton.Text = 'Plot';

            % Create PlotSettingPanel
            app.PlotSettingPanel = uipanel(app.PlotTab);
            app.PlotSettingPanel.TitlePosition = 'centertop';
            app.PlotSettingPanel.Title = 'Plot Setting';
            app.PlotSettingPanel.Position = [8 228 408 185];

            % Create FolderButton
            app.FolderButton = uibutton(app.PlotSettingPanel, 'push');
            app.FolderButton.ButtonPushedFcn = createCallbackFcn(app, @FolderButtonPushed, true);
            app.FolderButton.Position = [18 132 52 22];
            app.FolderButton.Text = 'Folder';

            % Create EditField
            app.EditField = uieditfield(app.PlotSettingPanel, 'text');
            app.EditField.Position = [90 132 251 22];

            % Create PlotfilesfromEditFieldLabel
            app.PlotfilesfromEditFieldLabel = uilabel(app.PlotSettingPanel);
            app.PlotfilesfromEditFieldLabel.HorizontalAlignment = 'right';
            app.PlotfilesfromEditFieldLabel.Position = [1 87 82 22];
            app.PlotfilesfromEditFieldLabel.Text = ' Plot files from';

            % Create PlotfilesfromEditField
            app.PlotfilesfromEditField = uieditfield(app.PlotSettingPanel, 'numeric');
            app.PlotfilesfromEditField.ValueDisplayFormat = '%.0f';
            app.PlotfilesfromEditField.Position = [90 87 48 22];

            % Create ToEditFieldLabel
            app.ToEditFieldLabel = uilabel(app.PlotSettingPanel);
            app.ToEditFieldLabel.HorizontalAlignment = 'right';
            app.ToEditFieldLabel.Position = [137 87 25 22];
            app.ToEditFieldLabel.Text = 'To';

            % Create ToEditField
            app.ToEditField = uieditfield(app.PlotSettingPanel, 'numeric');
            app.ToEditField.ValueDisplayFormat = '%.0f';
            app.ToEditField.Position = [172 87 48 22];

            % Create ResolutionDropDownLabel
            app.ResolutionDropDownLabel = uilabel(app.PlotSettingPanel);
            app.ResolutionDropDownLabel.HorizontalAlignment = 'right';
            app.ResolutionDropDownLabel.Position = [223 87 62 22];
            app.ResolutionDropDownLabel.Text = 'Resolution';

            % Create ResolutionDropDown
            app.ResolutionDropDown = uidropdown(app.PlotSettingPanel);
            app.ResolutionDropDown.Items = {'1', '2', ''};
            app.ResolutionDropDown.Position = [299 87 71 22];
            app.ResolutionDropDown.Value = '2';

            % Create NumberofpixelstoneglectattubesendEditFieldLabel
            app.NumberofpixelstoneglectattubesendEditFieldLabel = uilabel(app.PlotSettingPanel);
            app.NumberofpixelstoneglectattubesendEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofpixelstoneglectattubesendEditFieldLabel.Position = [15 17 228 22];
            app.NumberofpixelstoneglectattubesendEditFieldLabel.Text = 'Number of  pixels to neglect at tube''s end';

            % Create NumberofpixelstoneglectattubesendEditField
            app.NumberofpixelstoneglectattubesendEditField = uieditfield(app.PlotSettingPanel, 'numeric');
            app.NumberofpixelstoneglectattubesendEditField.ValueDisplayFormat = '%.0f';
            app.NumberofpixelstoneglectattubesendEditField.Position = [275 17 95 22];

            % Create PlotdfromEditFieldLabel
            app.PlotdfromEditFieldLabel = uilabel(app.PlotSettingPanel);
            app.PlotdfromEditFieldLabel.HorizontalAlignment = 'right';
            app.PlotdfromEditFieldLabel.Position = [-3 52 86 22];
            app.PlotdfromEditFieldLabel.Text = ' Plot d from (Å)';

            % Create PlotdfromEditField
            app.PlotdfromEditField = uieditfield(app.PlotSettingPanel, 'numeric');
            app.PlotdfromEditField.ValueDisplayFormat = '%.2f';
            app.PlotdfromEditField.Position = [90 52 48 22];

            % Create ToEditField_2Label
            app.ToEditField_2Label = uilabel(app.PlotSettingPanel);
            app.ToEditField_2Label.HorizontalAlignment = 'right';
            app.ToEditField_2Label.Position = [137 52 25 22];
            app.ToEditField_2Label.Text = 'To';

            % Create ToEditField_2
            app.ToEditField_2 = uieditfield(app.PlotSettingPanel, 'numeric');
            app.ToEditField_2.ValueDisplayFormat = '%.2f';
            app.ToEditField_2.Position = [172 52 48 22];

            % Create PoleFigureSettingPanel
            app.PoleFigureSettingPanel = uipanel(app.PlotTab);
            app.PoleFigureSettingPanel.TitlePosition = 'centertop';
            app.PoleFigureSettingPanel.Title = 'Pole Figure Setting';
            app.PoleFigureSettingPanel.Position = [9 34 410 155];

            % Create RotatearoundXdegreesEditFieldLabel
            app.RotatearoundXdegreesEditFieldLabel = uilabel(app.PoleFigureSettingPanel);
            app.RotatearoundXdegreesEditFieldLabel.HorizontalAlignment = 'right';
            app.RotatearoundXdegreesEditFieldLabel.Position = [61 90 151 22];
            app.RotatearoundXdegreesEditFieldLabel.Text = 'Rotate around X: (degrees)';

            % Create RotatearoundXdegreesEditField
            app.RotatearoundXdegreesEditField = uieditfield(app.PoleFigureSettingPanel, 'numeric');
            app.RotatearoundXdegreesEditField.ValueDisplayFormat = '%.1f';
            app.RotatearoundXdegreesEditField.Position = [227 90 100 22];

            % Create RotatearoundYdegreesEditFieldLabel
            app.RotatearoundYdegreesEditFieldLabel = uilabel(app.PoleFigureSettingPanel);
            app.RotatearoundYdegreesEditFieldLabel.HorizontalAlignment = 'right';
            app.RotatearoundYdegreesEditFieldLabel.Position = [62 47 150 22];
            app.RotatearoundYdegreesEditFieldLabel.Text = 'Rotate around Y: (degrees)';

            % Create RotatearoundYdegreesEditField
            app.RotatearoundYdegreesEditField = uieditfield(app.PoleFigureSettingPanel, 'numeric');
            app.RotatearoundYdegreesEditField.ValueDisplayFormat = '%.1f';
            app.RotatearoundYdegreesEditField.Position = [227 47 100 22];

            % Create GeneratePoleFigureButton
            app.GeneratePoleFigureButton = uibutton(app.PoleFigureSettingPanel, 'push');
            app.GeneratePoleFigureButton.ButtonPushedFcn = createCallbackFcn(app, @GeneratePoleFigureButtonPushed, true);
            app.GeneratePoleFigureButton.Position = [130 14 130 22];
            app.GeneratePoleFigureButton.Text = 'Generate Pole Figure';

            % Create AnalyseTab
            app.AnalyseTab = uitab(app.TabGroup);
            app.AnalyseTab.Title = 'Analyse';

            % Create UIAxes_3
            app.UIAxes_3 = uiaxes(app.AnalyseTab);
            xlabel(app.UIAxes_3, '\Phi')
            ylabel(app.UIAxes_3, {'\Theta'; ''})
            app.UIAxes_3.XTickLabelRotation = 0;
            app.UIAxes_3.YTickLabelRotation = 0;
            app.UIAxes_3.ZTickLabelRotation = 0;
            app.UIAxes_3.Position = [436 311 321 252];

            % Create AngleanalysisPanel
            app.AngleanalysisPanel = uipanel(app.AnalyseTab);
            app.AngleanalysisPanel.TitlePosition = 'centertop';
            app.AngleanalysisPanel.Title = 'Angle analysis';
            app.AngleanalysisPanel.Position = [15 326 393 230];

            % Create Thetaofpoint1EditFieldLabel
            app.Thetaofpoint1EditFieldLabel = uilabel(app.AngleanalysisPanel);
            app.Thetaofpoint1EditFieldLabel.HorizontalAlignment = 'right';
            app.Thetaofpoint1EditFieldLabel.Position = [13 116 89 22];
            app.Thetaofpoint1EditFieldLabel.Text = 'Theta of point 1';

            % Create Thetaofpoint1EditField
            app.Thetaofpoint1EditField = uieditfield(app.AngleanalysisPanel, 'numeric');
            app.Thetaofpoint1EditField.ValueDisplayFormat = '%.1f';
            app.Thetaofpoint1EditField.Position = [130 116 44 22];

            % Create Phiofpoint1EditFieldLabel
            app.Phiofpoint1EditFieldLabel = uilabel(app.AngleanalysisPanel);
            app.Phiofpoint1EditFieldLabel.HorizontalAlignment = 'right';
            app.Phiofpoint1EditFieldLabel.Position = [223 116 76 22];
            app.Phiofpoint1EditFieldLabel.Text = 'Phi of point 1';

            % Create Phiofpoint1EditField
            app.Phiofpoint1EditField = uieditfield(app.AngleanalysisPanel, 'numeric');
            app.Phiofpoint1EditField.ValueDisplayFormat = '%.1f';
            app.Phiofpoint1EditField.Position = [327 116 44 22];

            % Create Thetaofpoint2EditFieldLabel
            app.Thetaofpoint2EditFieldLabel = uilabel(app.AngleanalysisPanel);
            app.Thetaofpoint2EditFieldLabel.HorizontalAlignment = 'right';
            app.Thetaofpoint2EditFieldLabel.Position = [13 74 89 22];
            app.Thetaofpoint2EditFieldLabel.Text = 'Theta of point 2';

            % Create Thetaofpoint2EditField
            app.Thetaofpoint2EditField = uieditfield(app.AngleanalysisPanel, 'numeric');
            app.Thetaofpoint2EditField.ValueDisplayFormat = '%.1f';
            app.Thetaofpoint2EditField.Position = [130 74 44 22];

            % Create Phiofpoint2EditFieldLabel
            app.Phiofpoint2EditFieldLabel = uilabel(app.AngleanalysisPanel);
            app.Phiofpoint2EditFieldLabel.HorizontalAlignment = 'right';
            app.Phiofpoint2EditFieldLabel.Position = [223 74 76 22];
            app.Phiofpoint2EditFieldLabel.Text = 'Phi of point 2';

            % Create Phiofpoint2EditField
            app.Phiofpoint2EditField = uieditfield(app.AngleanalysisPanel, 'numeric');
            app.Phiofpoint2EditField.ValueDisplayFormat = '%.1f';
            app.Phiofpoint2EditField.Position = [327 74 44 22];

            % Create ReturnanglesButton
            app.ReturnanglesButton = uibutton(app.AngleanalysisPanel, 'push');
            app.ReturnanglesButton.ButtonPushedFcn = createCallbackFcn(app, @ReturnanglesButtonPushed, true);
            app.ReturnanglesButton.Position = [61 25 100 22];
            app.ReturnanglesButton.Text = 'Return angles';

            % Create EditField_3
            app.EditField_3 = uieditfield(app.AngleanalysisPanel, 'numeric');
            app.EditField_3.Position = [270 25 100 22];

            % Create LoadButton_2
            app.LoadButton_2 = uibutton(app.AngleanalysisPanel, 'push');
            app.LoadButton_2.ButtonPushedFcn = createCallbackFcn(app, @LoadButton_2Pushed, true);
            app.LoadButton_2.Position = [32 171 100 22];
            app.LoadButton_2.Text = 'Load';

            % Create FileinspectionPanel
            app.FileinspectionPanel = uipanel(app.AnalyseTab);
            app.FileinspectionPanel.TitlePosition = 'centertop';
            app.FileinspectionPanel.Title = 'File inspection';
            app.FileinspectionPanel.Position = [15 35 412 251];

            % Create ThetaminEditFieldLabel
            app.ThetaminEditFieldLabel = uilabel(app.FileinspectionPanel);
            app.ThetaminEditFieldLabel.HorizontalAlignment = 'right';
            app.ThetaminEditFieldLabel.Position = [28 175 59 22];
            app.ThetaminEditFieldLabel.Text = 'Theta min';

            % Create ThetaminEditField
            app.ThetaminEditField = uieditfield(app.FileinspectionPanel, 'numeric');
            app.ThetaminEditField.ValueDisplayFormat = '%.0f';
            app.ThetaminEditField.Position = [120 175 35 22];

            % Create ThetamaxEditFieldLabel
            app.ThetamaxEditFieldLabel = uilabel(app.FileinspectionPanel);
            app.ThetamaxEditFieldLabel.HorizontalAlignment = 'right';
            app.ThetamaxEditFieldLabel.Position = [198 175 62 22];
            app.ThetamaxEditFieldLabel.Text = 'Theta max';

            % Create ThetamaxEditField
            app.ThetamaxEditField = uieditfield(app.FileinspectionPanel, 'numeric');
            app.ThetamaxEditField.ValueDisplayFormat = '%.0f';
            app.ThetamaxEditField.Position = [293 175 35 22];

            % Create PhiminEditFieldLabel
            app.PhiminEditFieldLabel = uilabel(app.FileinspectionPanel);
            app.PhiminEditFieldLabel.HorizontalAlignment = 'right';
            app.PhiminEditFieldLabel.Position = [41 136 46 22];
            app.PhiminEditFieldLabel.Text = 'Phi min';

            % Create PhiminEditField
            app.PhiminEditField = uieditfield(app.FileinspectionPanel, 'numeric');
            app.PhiminEditField.ValueDisplayFormat = '%.0f';
            app.PhiminEditField.Position = [120 136 35 22];

            % Create PhimaxEditFieldLabel
            app.PhimaxEditFieldLabel = uilabel(app.FileinspectionPanel);
            app.PhimaxEditFieldLabel.HorizontalAlignment = 'right';
            app.PhimaxEditFieldLabel.Position = [211 136 49 22];
            app.PhimaxEditFieldLabel.Text = 'Phi max';

            % Create PhimaxEditField
            app.PhimaxEditField = uieditfield(app.FileinspectionPanel, 'numeric');
            app.PhimaxEditField.ValueDisplayFormat = '%.0f';
            app.PhimaxEditField.Position = [293 136 35 22];

            % Create GetfilesButton
            app.GetfilesButton = uibutton(app.FileinspectionPanel, 'push');
            app.GetfilesButton.ButtonPushedFcn = createCallbackFcn(app, @GetfilesButtonPushed, true);
            app.GetfilesButton.Position = [130 71 100 22];
            app.GetfilesButton.Text = 'Get files';

            % Create UITable
            app.UITable = uitable(app.AnalyseTab);
            app.UITable.ColumnName = {'File name'; 'max phi'; 'min phi'; 'max theta'; 'min theta'};
            app.UITable.RowName = {};
            app.UITable.FontSize = 14;
            app.UITable.Position = [456 35 301 250];

            % Create EditField_2
            app.EditField_2 = uieditfield(app.AnalyseTab, 'text');
            app.EditField_2.Position = [169 497 217 22];

            % Create AdvancedanalysisTab
            app.AdvancedanalysisTab = uitab(app.TabGroup);
            app.AdvancedanalysisTab.Title = 'Advanced analysis';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.AdvancedanalysisTab);
            xlabel(app.UIAxes2, '\Phi')
            ylabel(app.UIAxes2, '\Theta')
            app.UIAxes2.XTickLabelRotation = 0;
            app.UIAxes2.YTickLabelRotation = 0;
            app.UIAxes2.ZTickLabelRotation = 0;
            app.UIAxes2.Position = [389 343 353 243];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.AdvancedanalysisTab);
            title(app.UIAxes3, 'Mosaic Angle')
            xlabel(app.UIAxes3, '\Delta{\Theta}')
            ylabel(app.UIAxes3, 'Intensity')
            app.UIAxes3.XTickLabelRotation = 0;
            app.UIAxes3.YTickLabelRotation = 0;
            app.UIAxes3.ZTickLabelRotation = 0;
            app.UIAxes3.Position = [389 63 353 257];

            % Create LoadandPlotPanel
            app.LoadandPlotPanel = uipanel(app.AdvancedanalysisTab);
            app.LoadandPlotPanel.TitlePosition = 'centertop';
            app.LoadandPlotPanel.Title = 'Load and Plot';
            app.LoadandPlotPanel.Position = [12 354 368 220];

            % Create LoadButton
            app.LoadButton = uibutton(app.LoadandPlotPanel, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Position = [12 154 100 22];
            app.LoadButton.Text = 'Load';

            % Create EditField_4
            app.EditField_4 = uieditfield(app.LoadandPlotPanel, 'text');
            app.EditField_4.Position = [152 154 157 22];

            % Create IntensitycutonDropDownLabel
            app.IntensitycutonDropDownLabel = uilabel(app.LoadandPlotPanel);
            app.IntensitycutonDropDownLabel.HorizontalAlignment = 'right';
            app.IntensitycutonDropDownLabel.Position = [14 116 86 22];
            app.IntensitycutonDropDownLabel.Text = 'Intensity cut on';

            % Create IntensitycutonDropDown
            app.IntensitycutonDropDown = uidropdown(app.LoadandPlotPanel);
            app.IntensitycutonDropDown.Items = {'Theta', 'Phi'};
            app.IntensitycutonDropDown.Position = [154 116 157 22];
            app.IntensitycutonDropDown.Value = 'Theta';

            % Create ThetaSliderLabel
            app.ThetaSliderLabel = uilabel(app.LoadandPlotPanel);
            app.ThetaSliderLabel.HorizontalAlignment = 'right';
            app.ThetaSliderLabel.Position = [12 35 36 22];
            app.ThetaSliderLabel.Text = 'Theta';

            % Create ThetaSlider
            app.ThetaSlider = uislider(app.LoadandPlotPanel);
            app.ThetaSlider.Limits = [0 180];
            app.ThetaSlider.ValueChangingFcn = createCallbackFcn(app, @ThetaSliderValueChanging, true);
            app.ThetaSlider.Position = [69 44 173 3];

            % Create ThetaEditField
            app.ThetaEditField = uieditfield(app.LoadandPlotPanel, 'numeric');
            app.ThetaEditField.Position = [272 25 41 22];

            % Create PhiEditField
            app.PhiEditField = uieditfield(app.LoadandPlotPanel, 'numeric');
            app.PhiEditField.Position = [271 73 41 22];

            % Create PhiLabel
            app.PhiLabel = uilabel(app.LoadandPlotPanel);
            app.PhiLabel.HorizontalAlignment = 'right';
            app.PhiLabel.Position = [14 85 25 22];
            app.PhiLabel.Text = 'Phi';

            % Create PhiSlider
            app.PhiSlider = uislider(app.LoadandPlotPanel);
            app.PhiSlider.Limits = [0 355];
            app.PhiSlider.ValueChangingFcn = createCallbackFcn(app, @PhiSliderValueChanging, true);
            app.PhiSlider.Position = [60 94 184 3];

            % Create GrainclassificationPanel
            app.GrainclassificationPanel = uipanel(app.AdvancedanalysisTab);
            app.GrainclassificationPanel.TitlePosition = 'centertop';
            app.GrainclassificationPanel.Title = 'Grain classification';
            app.GrainclassificationPanel.Position = [15 30 365 314];

            % Create ThresholdEditFieldLabel
            app.ThresholdEditFieldLabel = uilabel(app.GrainclassificationPanel);
            app.ThresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.ThresholdEditFieldLabel.Position = [11 252 59 22];
            app.ThresholdEditFieldLabel.Text = 'Threshold';

            % Create ThresholdEditField
            app.ThresholdEditField = uieditfield(app.GrainclassificationPanel, 'numeric');
            app.ThresholdEditField.Position = [89 252 46 22];

            % Create MosaicrangeEditFieldLabel
            app.MosaicrangeEditFieldLabel = uilabel(app.GrainclassificationPanel);
            app.MosaicrangeEditFieldLabel.HorizontalAlignment = 'right';
            app.MosaicrangeEditFieldLabel.Position = [163 252 78 22];
            app.MosaicrangeEditFieldLabel.Text = 'Mosaic range';

            % Create MosaicrangeEditField
            app.MosaicrangeEditField = uieditfield(app.GrainclassificationPanel, 'numeric');
            app.MosaicrangeEditField.Position = [260 252 46 22];

            % Create DetectpeaksButton
            app.DetectpeaksButton = uibutton(app.GrainclassificationPanel, 'push');
            app.DetectpeaksButton.ButtonPushedFcn = createCallbackFcn(app, @DetectpeaksButtonPushed, true);
            app.DetectpeaksButton.Position = [131 221 100 22];
            app.DetectpeaksButton.Text = 'Detect peaks';

            % Create UITable2
            app.UITable2 = uitable(app.GrainclassificationPanel);
            app.UITable2.ColumnName = {'Row number'; 'Theta'; 'Phi'; 'Counts'};
            app.UITable2.RowName = {};
            app.UITable2.Position = [22 93 317 121];

            % Create RownumberEditFieldLabel
            app.RownumberEditFieldLabel = uilabel(app.GrainclassificationPanel);
            app.RownumberEditFieldLabel.HorizontalAlignment = 'right';
            app.RownumberEditFieldLabel.Position = [35 60 74 22];
            app.RownumberEditFieldLabel.Text = 'Row number';

            % Create RownumberEditField
            app.RownumberEditField = uieditfield(app.GrainclassificationPanel, 'numeric');
            app.RownumberEditField.ValueDisplayFormat = '%.0f';
            app.RownumberEditField.Position = [163 60 43 22];

            % Create DeleteButton
            app.DeleteButton = uibutton(app.GrainclassificationPanel, 'push');
            app.DeleteButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteButtonPushed, true);
            app.DeleteButton.Position = [249 60 100 22];
            app.DeleteButton.Text = 'Delete';

            % Create TheoreticalanglesEditFieldLabel
            app.TheoreticalanglesEditFieldLabel = uilabel(app.GrainclassificationPanel);
            app.TheoreticalanglesEditFieldLabel.HorizontalAlignment = 'right';
            app.TheoreticalanglesEditFieldLabel.Position = [35 33 104 22];
            app.TheoreticalanglesEditFieldLabel.Text = 'Theoretical angles';

            % Create TheoreticalanglesEditField
            app.TheoreticalanglesEditField = uieditfield(app.GrainclassificationPanel, 'numeric');
            app.TheoreticalanglesEditField.Position = [163 33 43 22];

            % Create GrainanalysisButton
            app.GrainanalysisButton = uibutton(app.GrainclassificationPanel, 'push');
            app.GrainanalysisButton.ButtonPushedFcn = createCallbackFcn(app, @GrainanalysisButtonPushed, true);
            app.GrainanalysisButton.Position = [249 33 100 22];
            app.GrainanalysisButton.Text = 'Grain analysis';

            % Create BasisnumberEditFieldLabel
            app.BasisnumberEditFieldLabel = uilabel(app.GrainclassificationPanel);
            app.BasisnumberEditFieldLabel.HorizontalAlignment = 'right';
            app.BasisnumberEditFieldLabel.Position = [35 1 79 22];
            app.BasisnumberEditFieldLabel.Text = 'Basis number';

            % Create BasisnumberEditField
            app.BasisnumberEditField = uieditfield(app.GrainclassificationPanel, 'numeric');
            app.BasisnumberEditField.ValueDisplayFormat = '%.0f';
            app.BasisnumberEditField.Position = [134 2 43 22];

            % Create EditField_5
            app.EditField_5 = uieditfield(app.GrainclassificationPanel, 'numeric');
            app.EditField_5.ValueDisplayFormat = '%.0f';
            app.EditField_5.Position = [205 2 40 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GUI6_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
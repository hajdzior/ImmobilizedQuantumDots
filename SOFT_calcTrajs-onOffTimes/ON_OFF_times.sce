// ON-OFF analysis script

//------------------------------
// Last mod: 17-07-14

// Calculates:
// - the thresholds
// - on and off times

// input folder: 'trajsTXT'
// output folder: 'on_off_times'
// MW 17-07-14
//-------------------------------

//#################################
//                                #
//            FUNCTIONS           #
//                                #
//#################################
// 
// //fit_on-function
//  function [e]=G(p,z)
//    if z(1) > 0 then
//      e = (z(2)-size_on_off(1)*h*exp(-z(1)/p(1))/p(1))./sqrt(z(1)) 
//    else
//      e = 0
//    end //if
//  endfunction

// // LS fit function 
// function y = FF(hist_i, p)
//     y=sample_size*exp(-p(1).*(edges_left(hist_i))).*(1.0-exp(-p(1).*sizes(hist_i)))
// endfunction

// // derivative
// function y = dFF(hist_i, p)
//     y=sample_size*exp(-p(1).*(edges_left(hist_i))).*(-(edges_left(hist_i))+(edges_right(hist_i)).*exp(-p(1)*sizes(hist_i)))
// endfunction

stacksize(50000000)

// LS weigthed cost function 
function e = LSW(p, z),
    y = z(1), x = z(2);
    if ( y > hist_cutoff ) then 
      e = (y - FF(x, p)) / sqrt(y)
      else e = 0.0;
    end
endfunction

// histogram fit
function [ k_fit, sigma_fit, chisq_fit ] = histfit(LSW, FF, dFF, Z, k_ini, edges, hist_cutoff)
    edges_left = edges( 1:$-1 );
    edges_right = edges( 2:$ );
    sizes = edges_right - edges_left;
    sample_size = sum( Z(1,:) );
    [ k_fit, err_fit ] = datafit( LSW, Z, k_ini );
    chisq_fit = err_fit / ( length( find( Z(1,:) > hist_cutoff ) ) -1 );
    
    F = zeros( length( Z(1,:) ), 1);
    for i = 1:length( Z(1,:) ), 
        if ( Z(1,i) > hist_cutoff ) F(i,1) = dFF(i, k_fit) ./ sqrt( Z(1,i) );
        else F(i,1) = 0.0;
        end
    end  

    cov_fit = inv( F'*F )
    sigma_fit = sqrt( cov_fit )
endfunction

//thresholds
 function [cutoff, cutoff_high_3, cutoff_low_3, cutoff_high_5, cutoff_low_5] = calc_cutoffs( events )
   // initial estimate
   cutoff = ( max(events) + min(events) ) / 2.0; 
   while(1)
      avg_off = sum( events(find(events<cutoff)) ) / length(events(find(events<cutoff)));
      avg_on = sum( events(find(events>cutoff)) ) / length(events(find(events>cutoff)));
      cutoff_new = ( avg_off + avg_on ) / 2.0
      if abs( cutoff_new - cutoff ) / cutoff < 0.00001 then break;
      else cutoff = cutoff_new;
      end
    end;//while
    cutoff_low_3 = avg_off + (1/3.0)*(avg_on - avg_off);
    cutoff_high_3 = avg_off + (2/3.0)*(avg_on - avg_off);
    cutoff_low_5 = avg_off + (1/5.0)*(avg_on - avg_off);
    cutoff_high_5 = avg_off + (4/5.0)*(avg_on - avg_off);
 endfunction

// derivative $$$ need to standard dev
function y = dFF(hist_i,p)
    y=sample_size*exp(-p(1).*(edges_left(hist_i))).*(-(edges_left(hist_i))+(edges_right(hist_i)).*exp(-p(1)*sizes(hist_i)))
endfunction

//vizualization, theoretical curve and fit
function  [err, occ, x_mid, h, n]=histogram(on_off, binL, n)
 //"x" for histograms
 maxim = max(on_off);
 minim = min(on_off);
 size_on_off = size(on_off);
//bin size for histogram
 //n = 15; //start numbers of bins
 h_temp = (maxim - minim)/n
 h = (ceil(h_temp/binL)+1)*binL;
 x_mid=[minim+h/2:h:n*h+h/2]';
 val = linspace(minim,n*h,n+1)'
 [err, occ] = dsearch(on_off, val)
endfunction

// //cumulative histograms
// function  [cum_hist]=cumulative_histogram(on_off_occ)
//   tmp_occ = 0;
//   occ_size = size(on_off_occ)
//   for(ch = 1:1:occ_size(1))
//     occ_sum = sum(on_off_occ)
//     tmp_occ = tmp_occ + on_off_occ(ch)
//     cum_hist(ch) = tmp_occ/occ_sum
//   end  
// endfunction


PLOTS_FACTOR = 0;              // make plots

// SIMULATION = 0;                // run simulation
// 
// PCH_analysis         = 0;       
// ICF_analysis         = 0;
ON_OFF_analysis      = 1;
//In On-OFF analysis:
trajectory       = 1;
  
thresh_1             = 1;
thresh_2_1_3        = 1;
thresh_2_1_5        = 1; 
thresh_2_st_dev     = 0;
thresh_custom       = 0;

//#################################
//                                #
//           Parameters           #
//      needed by the fitting     #
//           algorithm            #
//                                #
//#################################   

// k_on_ini = 0.5
// k_off_ini = 0.5
// I_on = 40
// I_off = 10
thresh_custom_up = 30
thresh_custom_down = 20

//#################################
//                                #
//   Now everything is ready      #
//   program can start analysis   #
//                                #
//#################################   

WHICH_DATA = [trajectory];
INPUT_DATA = ['../tmp/trajectory.txt'];
INPUT_DATA_for_filenames = ['czyDziala']
my_thresholds = ['thresh-1', 'thresh-2--1-3', 'thresh-2--1-5', 'thresh_2_st_dev', 'thresh_custom'];
onoffANALYSIS = [thresh_1, thresh_2_1_3, thresh_2_1_5, thresh_2_st_dev, thresh_custom];

make_analysis = 0;

data_in = INPUT_DATA;
data = fscanfMat(data_in);
time = data(:,1); events = data(:,2);
binL = time(2)-time(1)
thr = mopen('../tmp/thresholds.txt', 'w')
thr_v = mopen('../tmp/thresholds_val.txt', 'w')
thr_n = mopen('../tmp/thresholds_names.txt', 'w')
//thresholds
[cutoff, cutoff_high_3, cutoff_low_3, cutoff_high_5, cutoff_low_5] = calc_cutoffs( events );
for(thresh = 1:1:5)  
	if (thresh == 1 & onoffANALYSIS(thresh) == 1)
		thresh_u = cutoff;
		thresh_d = cutoff;
		make_analysis = 1;
	end

	if (thresh == 2 & onoffANALYSIS(thresh) == 1) 
		thresh_u = cutoff_high_3;
		thresh_d = cutoff_low_3;
		make_analysis = 1;
	end

	if (thresh == 3 & onoffANALYSIS(thresh) == 1)
		thresh_u = cutoff_high_5
		thresh_d = cutoff_low_5
		make_analysis = 1;
	end

	if (thresh == 4 & onoffANALYSIS(thresh) == 1)
		if (I_on > I_off)
			thresh_u = I_on - I_on^(0.5)
			thresh_d = I_off + I_off^(0.5)
		elseif(I_off > I_on)
			thresh_u = I_off - I_off^(0.5)
			thresh_d = I_on + I_on^(0.5)
		end
		make_analysis = 1;
	end

	if (thresh == 5 & onoffANALYSIS(thresh) == 1)
		thresh_u = thresh_custom_up
		thresh_d = thresh_custom_down
		make_analysis = 1;
	end

	if (make_analysis == 1)
		events_size = size(events);
		states = [1:events_size(1)];

		//initial state
		if (events(1) > cutoff) then state = 1;
		elseif (events(1) <= cutoff) then state = 2;
		end //if

		//"on" and "off" states
		for i = 1:events_size(1)//*
			if (events(i) > thresh_u) then state = 1;
			elseif (events(i) < thresh_d) then state = 2;
			end //if
			states(i) = state;
		end //for*

		//making states data
		threshold_now = my_thresholds(thresh);
		input_data_now =  INPUT_DATA_for_filenames
		my_times_on = 'fit_results/ON_times_' +  '_' + threshold_now + '.txt';
		my_times_off = 'fit_results/OFF_times_' + '_' + threshold_now + '.txt';
		states_size = size(states);
		volume = 1;
		k=1;
		l=1;
		on = zeros();
		off = zeros();
		for j = 1:events_size(1)-1
			if(states(j) == states(j+1)) then volume = volume + 1;
			else time_v = volume * binL;
				if(states(j) == 1) then on(k) = time_v;, k = k + 1;
				elseif(states(j) == 2) then  off(l) =  time_v;, l = l+1;
				end
				volume = 1;
			end
		end

// 		
// 
// 	//#################################
// 	//                                #
// 	//              Fit               #
// 	//                                #
// 	//#################################
// 
// 		n = 15 //number of histograms bins
// 		//ON
// 		which = 1;
// 		[on_err, on_occ, on_x_mid, on_h]=histogram(on, binL, n)
// 		Z_on=[on_occ';1:length(on_occ)'];
// 		range_on = n;
// 		edges_on = (0:range_on) * on_h;
// 		hist_cuton = 0;
// 		//fit
// 		[ k_on_fit, sigma_on_fit, chisq_on_fit] = histfit( LSW, FF, dFF, Z_on, k_on_ini, edges_on, hist_cuton );
// 		tau_fit = 1/k_on_fit;
// 		sample_size = sum( on_occ );
// 		edges_left = edges_on( 1:$-1 );
// 		edges_right = edges_on( 2:$ );
// 		sizes = edges_right - edges_left;
// 
// 	
// 		//OFF
// 		which = 2;
// 		[off_err, off_occ, off_x_mid, off_h]=histogram(off, binL, n);
// 		Z_off=[off_occ';1:length(off_occ)'];
// 		range_off = n;
// 		edges_off = (0:range_off) * off_h;
// 		hist_cutoff = 0;
// 		//fit
// 		[ k_off_fit, sigma_off_fit, chisq_off_fit] = histfit( LSW, FF, dFF, Z_off, k_off_ini, edges_off, hist_cutoff );
// 		tau_fit = 1/k_off_fit;
// 		sample_size = sum( off_occ );
// 		edges_left = edges_off( 1:$-1 );
// 		edges_right = edges_off( 2:$ );
// 		sizes = edges_right - edges_left;
// 
// 		//writing to the file
		OnOff_times_f = '../tmp/on-off_times'  + '_' + threshold_now + '.txt'
// 		off_times_f = '../tmp/off_times'  + '_' + threshold_now + '.txt'
// 		thresh_f = 'thresh/thr'  + '_' + threshold_now + '.txt'
		fd_OnOff_t = mopen(OnOff_times_f, 'w')
// 		fd_off_t = mopen(off_times_f, 'w')
		thr_fil = threshold_now + ' & %16.2f \t &%16.2f \\\\ \n';
		mfprintf(thr, thr_fil, thresh_u, thresh_d)
		mfprintf(thr_n, '%s \n', threshold_now)
		mfprintf(fd_OnOff_t, '%16.8f \t %16.8f \n', on(:,:), off(:,:))
// 		mfprintf(fd_off_t, '%16.8f \n', off(:,:))
		mclose(fd_OnOff_t);
// 		mclose(fd_off_t);

// 		//writing to the file
		OnOff_times_fVal = '../tmp/on-off_times'  + '_val_' + threshold_now + '.txt'
// 		off_times_f = '../tmp/off_times'  + '_' + threshold_now + '.txt'
// 		thresh_f = 'thresh/thr'  + '_' + threshold_now + '.txt'
// 		fd_OnOff_tVal = mopen(OnOff_times_f, 'w')
// 		fd_off_t = mopen(off_times_f, 'w')
// 		thr_fil = threshold_now + ' & %16.2f \t &%16.2f \\\\ \n';
		mfprintf(thr_v, '%5.2f \t %5.2f \n', thresh_u, thresh_d)
// 		mfprintf(fd_off_t, '%16.8f \n', off(:,:))
// 		mclose(fd_off_t);

// 	//#################################
// 	//                                #
// 	//           Plots                #
// 	//                                #
// 	//#################################
// 			
// 		if PLOTS_FACTOR == 1
// 			//zmienic trzy kolejne linie (sa 2 razy)
// 			edges_left = edges_on( 1:$-1 );
// 			edges_right = edges_on( 2:$ );
// 			sizes = edges_right - edges_left;
// 			scf(3); clf(3);
// 			xtitle ("ON times fit")
// 			plot2d(on_x_mid, on_occ, style = -2)
// 			hist_on_fit = FF( 1:length(on_occ), k_on_fit );//fit curve
// 			plot2d( on_x_mid, hist_on_fit, 5 ) 
// 
// 
// 			scf(4); clf(4);
// 			edges_left = edges_off( 1:$-1 );
// 			edges_right = edges_off( 2:$ );
// 			sizes = edges_right - edges_left;
// 			xtitle ("OFF times fit")
// 			plot2d(off_x_mid, off_occ, style = -2)
// 			hist_off_fit = FF( 1:length(off_occ), k_off_fit );
// 			plot2d( off_x_mid, hist_off_fit, 5 )
//              
// 			//saving plots
// 			filename = 'fit_results/ON_analysis_' + input_data_now + '_' + threshold_now;
// 			xs2png(3,filename)
// 			filename = 'fit_results/OFF_analysis_' + input_data_now + '_' + threshold_now;
// 			xs2png(4,filename)
// 			xtitle ("OFF times cumulative histogram")
// 		end
// 			
 			make_analysis = 0; 
	end//if (make_analysis == 1)
end//for(thresh = 1:1:4)
// onoffANALYSIS(4) = 1;
//         end // if (WHICH_DATA(i_d) == 1)
//     end //for(i_d = 1:1:3)
mclose(thr);
mclose(thr_v);
mclose(thr_n);
exit;


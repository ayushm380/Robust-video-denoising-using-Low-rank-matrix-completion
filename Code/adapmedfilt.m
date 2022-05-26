function [denoised]=adapmedfilt(noisy, win_size)
    % Inputs:
    %   noisy : noisy image
    %   win_size  : maximum size of window
    % Outputs:
    %   denoised : denoised image


    denoised = noisy;
    process_done = false(size(noisy));

    % Iteration on k
    for k = 3:2:win_size
        
        min_pix_val = ordfilt2(noisy, 1, ones(k, k), 'symmetric');              % Obtaining min_pix_val of every pixel within the window of size k
        max_pix_val = ordfilt2(noisy, k * k, ones(k, k), 'symmetric');          % Obtaining max_pix_val of every pixel within the window of size k
        med_pix_val = medfilt2(noisy, [k k], 'symmetric');                      % Obtaining med_pix_val of every pixel within the window of size k

        step_B_vals = (med_pix_val > min_pix_val) & (max_pix_val > med_pix_val) & ~process_done;       % pixels that will go to step B
        B_pix_val = (noisy > min_pix_val) & (max_pix_val > noisy);

        op = step_B_vals & B_pix_val;                                  % pixel position wrt true op of step AB
        op_med_val = step_B_vals & ~B_pix_val;                         % pixel position wrt true median of op that satisfies A but not B

        denoised(op) = noisy(op);                                      % Updating the pixel values
        denoised(op_med_val) = med_pix_val(op_med_val);

        process_done = process_done | step_B_vals;                     % Updating the subset of processed values

        if all(process_done(:))
            break;
        end
    end
    denoised(~process_done) = med_pix_val(~process_done);
end
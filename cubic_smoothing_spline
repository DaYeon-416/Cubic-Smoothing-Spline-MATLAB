%% Read Matrix
csv_name = 'Wage.csv';
matrix   = readmatrix(csv_name);

% grouping with age - wage
grouping_with_age = groupsummary(matrix, matrix(:, 3), 'median');
age_wage_array    = num2cell([grouping_with_age(:, 3), grouping_with_age(:, 13)]);

% data - Column Vector 사용
xdata = cell2mat(age_wage_array(:, 1)); % csv - age
ydata = cell2mat(age_wage_array(:, 2)); % csv - wage

% % some example data
% xdata = [0,1,2,4]';
% ydata = [0,1,0,1]';

%% Parameter Selection
weights = ones(size(xdata)) * 1; 
% weights(25:50) = 0.001;

p = 0.001; 

normalizedsmooth = false; 

%% Main & Plot Graph
[x, y, w]      = prepare_data(xdata, ydata, weights); % 입력 데이터 가공
[coef_list, p] = make_spline(x, y, w, p, normalizedsmooth); % smoothing parameter 및 3차 다항식 coefficient 계산

% -------------------- # Select Derivative # --------------------
%  - 0 : no derivative
%  - 1 : first derivative
%  - 2 : second derivative
% ---------------------------------------------------------------
nu = 0;

break_point_list = x;
x_list = linspace(min(xdata), max(xdata), 70);
y_list = construct_fast(coef_list, break_point_list, x_list, nu); % x값과 coefficient를 통해 y값 계산

% plot
figure('Renderer', 'painters', 'Position', [800 200 900 600])
scatter(xdata, ydata, 50, 'MarkerFaceColor', "blue", 'MarkerEdgeColor', [0.3010 0.7450 0.9330]); hold on; grid on;
xlabel('x'); ylabel('y'); title('Cubic Smoothing Spline')
plot(x_list, y_list, 'o-', 'MarkerSize', 7, 'LineWidth', 1.2)

%% Function - prepare_data
function [x, y, w] = prepare_data(xdata, ydata, weights)
    x = (xdata)';
    y = (ydata)';
    
    if isnan(weights)
        w = ones(1, length(x));

    elseif any(weights(:) == 0)
        error_msg = 'Unable to use "0" in weight.';
        error(error_msg)

    else
        w = weights';
        if size(weights) ~= size(xdata)
            error('Weights size must be equal of xdata size')
        end
    end
end

%% Function - make_spline
function [coef_list, p] = make_spline(x, y, w, p, normalizedsmooth)

    if normalizedsmooth == false && ~isnan(p)
        fprintf('\nSmoothing Parameter "p" : %3f\n', p)
    end

    N  = size(x, 2);
    dx = diff(x);
    
    % Create diagonal sparse matrices
    % spdiags 함수로 만들어진 변수는 full()로 확인 가능
    diags_R  = vertcat(dx(2:end), 2*(dx(2:end) + dx(1:end-1)), dx(1:end-1));
    dx_denom = 1./dx;
    diags_QT = vertcat(dx_denom(1:end-1), -(dx_denom(2:end) + dx_denom(1:end-1)), dx_denom(2:end));
    
    R     = spdiags(diags_R', -1:1, N-2, N-2);
    QT    = spdiags(diags_QT', [0, 1, 2], N-2, N);
    D2    = spdiags((1./w)', 0, N, N);
    QTD2Q = QT * D2 * QT';

    if normalizedsmooth == true
       p = normalize_smooth(x, w, p);
       fprintf('\nNormalized Smoothing Parameter "p" : %3f\n', p)
    end
    
    if isnan(p)
       p = compute_smooth(R, QTD2Q);
       fprintf('\nAuto Smoothing Parameter "p" : %3f\n', p)
    end
    
    % Solve linear system for the 2nd derivatives
    a = 6.*(1.-p) * QTD2Q + p*R;
    b = QT * y';
    
    % Solve the sparse linear system ax=b
    u = a \ b; 
    dx = dx';
    
    D2Qu = D2 * QT' * u;
    yi = y' - 6.*(1.-p) * D2Qu;
    
    c = 3 * p * u;
    c_pad = zeros(size(c, 1) + 2, 1); % zero padding
    c_pad(2:end-1) = c;
    c = c_pad;

    c4 = yi(1:end-1, :); % a
    c3 = diff(yi)./dx - dx .* (2/3 * c(1:end-1, :) + 1/3 * c(2:end, :)); % b
    c2 = c(1:end-1, :); % c
    c1 = diff(c) ./ (3 * dx); % d
    
    c_shape = [N-1, 4];
    vert_c = vertcat(c1, c2, c3, c4);
    c = reshape(vert_c, c_shape)';
    
    coef_list = c;
end

%% Function - normalize_smooth
function p = normalize_smooth(x, w, p)

    span = max(x) - min(x);
    
    eff_x = 1 + (span.^2) / sum(diff(x).^2);
    eff_w = sum(w).^2 / sum(w.^2);
    k     = 80 * (span.^3) * (length(x).^(-2)) * (eff_x.^(-0.5)) * (eff_w.^(-0.5));
    
    if isnan(p)
        s = 0.5;
    else 
        s = p;
    end
    
    p = s/(s+(1-s)*k);
end

%% function - compute_smooth
function p = compute_smooth(R, QTD2Q)
    % -----------------------------------------------------------
    % (6*(1-p)*Q^T*D^2*Q + p*R)c  <-->  (6(1-p)*QTD2Q + p*R)c

    % (6(1-p)*QTD2Q + p*R) is a matrix of coefficient.
    % The Coefficient Matrix has a form of p*A + (1-p)*B
    % In case p is not provided, suppose equal weights are provided in matrix A and B.

    % So, in order to find p,
    % trace(6*(1-p)*QTD2Q) = trace(p*R)

    % if r = 6*trace(QTD2Q) / trace(R),
    % p / (1 - p) = r, so 
    % p = r / (1 + r)
    % -----------------------------------------------------------
    r = 6*trace(QTD2Q) / trace(R);
    p = r / (1 + r);
end

%% Function - construct_fast
function y_list = construct_fast(coef_list, break_point_list, x_list, nu)

    if length(break_point_list) < 2
        disp("error");
    end

    y_list = zeros(size(x_list));

    for i_x = 1:numel(x_list)
        x = x_list(i_x);
        canInterpolate = false;

        for i_b = 1:(length(break_point_list) - 1)

            if x >= break_point_list(i_b) && x < break_point_list(i_b + 1)
                dx = x - break_point_list(i_b);
                % y_list(i_x) = c4 + c3(x - xi) + c2(x - xi)^2 + c1(x - xi)^3
                % --> y_list(i_x) = c4 + dx*(c3 + dx*(c2 + dx*(c1)));
                % c1 : 3차항 계수 / c2 : 2차항 계수 / c3 : 1차항 계수 / c4 : 상수
                
                if nu == 0
                    y_list(i_x) = coef_list(4, i_b) + dx * (coef_list(3, i_b) + dx * (coef_list(2, i_b) + dx * (coef_list(1, i_b))));
                    break;
                elseif nu == 1
                    y_list(i_x) = coef_list(3, i_b) + dx * (2 * coef_list(2, i_b) + 3 * dx * coef_list(1, i_b));
                    break;
                elseif nu == 2
                    y_list(i_x) = 2 * coef_list(2, i_b) + 6 * dx * coef_list(1, i_b);
                end
            end
        end

        if i_b == size(coef_list, 2) && canInterpolate == false
            dx = x - break_point_list(i_b);
            
            if nu == 0
                y_list(i_x) = coef_list(4, i_b) + dx * (coef_list(3, i_b) + dx * (coef_list(2, i_b) + dx * (coef_list(1, i_b))));
            elseif nu == 1
                y_list(i_x) = coef_list(3, i_b) + dx * (2 * coef_list(2, i_b) + 3 * coef_list(1, i_b) * dx);
            elseif nu ==2
                y_list(i_x) = 2 * coef_list(2, i_b) + 6 * dx * coef_list(1, i_b);
            end
        end
    end
end

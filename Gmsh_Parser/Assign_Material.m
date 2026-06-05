function Mesh = Assign_Material(Mesh, omega, eps_0,omega_tilde,PML_set)
    % ---------------------------------------------------------
    % Assign material properties to each element
    %
    % Conductive materials are modeled using effective complex
    % relative permittivity:
    %
    % eps_r_eff = eps_r - 1j * sigma / (omega * eps_0)
    %
    % This convention is consistent with exp(j*w*t).
    % ---------------------------------------------------------

    if omega == 0
        error('omega must be nonzero for frequency-domain material assignment.');
    end

    % ---------------------------------------------------------
    % 1. Tensor helper functions
    % ---------------------------------------------------------
    iso_tensor = @(val) [val, 0, 0, ...
                         0, val, 0, ...
                         0, 0, val];

    % Conductive material tensor
    cond_eps_tensor = @(eps_r, sigma) ...
        iso_tensor(eps_r - 1j * sigma / (omega * eps_0));

    % ---------------------------------------------------------
    % 2. Material library
    % eps_r : relative permittivity
    % mu_r  : relative permeability
    % sigma : conductivity [S/m]
    % ---------------------------------------------------------
    MatLib = struct();

    % [1] vaccum and air
    MatLib.vaccum = struct( ...
        'eps_r', iso_tensor(1.0), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    MatLib.air = struct( ...
        'eps_r', iso_tensor(1.0006), ...
        'mu_r',  iso_tensor(1.0), ... 
        'sigma', 0.0);

    % [2] Dielectrics and polymers
    MatLib.teflon = struct( ...
        'eps_r', iso_tensor(2.1), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    MatLib.polyethylene = struct( ...
        'eps_r', iso_tensor(2.25), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    MatLib.nylon = struct( ...
        'eps_r', iso_tensor(4.0), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    MatLib.rubber = struct( ...
        'eps_r', iso_tensor(3.0), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    MatLib.wood_dry = struct( ...
        'eps_r', iso_tensor(2.0), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    % [3] PCB and RF substrates
    MatLib.rogers4003c = struct( ...
        'eps_r', iso_tensor(3.38), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    MatLib.fr4 = struct( ...
        'eps_r', iso_tensor(4.4), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    MatLib.alumina = struct( ...
        'eps_r', iso_tensor(9.4), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    % [4] Glass and ceramics
    MatLib.glass = struct( ...
        'eps_r', iso_tensor(4.2), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    MatLib.mica = struct( ...
        'eps_r', iso_tensor(6.0), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    % [5] Semiconductors and water
    MatLib.silicon = struct( ...
        'eps_r', iso_tensor(11.9), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    MatLib.gaas = struct( ...
        'eps_r', iso_tensor(12.9), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    MatLib.water = struct( ...
        'eps_r', iso_tensor(80.1), ...
        'mu_r',  iso_tensor(1.0), ...
        'sigma', 0.0);

    % ---------------------------------------------------------
    % [6] Metals
    %
    % At low frequency / RF / microwave range, metals are modeled
    % using conductivity sigma. The effective complex permittivity
    % is assigned to eps_r.
    % ---------------------------------------------------------
    MatLib.silver = struct( ...
        'eps_r', cond_eps_tensor(1.0, 6.30e7), ...
        'mu_r',  iso_tensor(0.99998), ...
        'sigma', 6.30e7);

    MatLib.gold = struct( ...
        'eps_r', cond_eps_tensor(1.0, 4.10e7), ...
        'mu_r',  iso_tensor(0.99996), ...
        'sigma', 4.10e7);

    MatLib.copper = struct( ...
        'eps_r', cond_eps_tensor(1.0, 5.96e7), ...
        'mu_r',  iso_tensor(0.99999), ...
        'sigma', 5.96e7);

    MatLib.aluminum = struct( ...
        'eps_r', cond_eps_tensor(1.0, 3.50e7), ...
        'mu_r',  iso_tensor(1.00002), ...
        'sigma', 3.50e7);
    MatLib.tungsten = struct( ...
        'eps_r', cond_eps_tensor(1.0, 1.79e7), ...
        'mu_r',  iso_tensor(1.00007), ...
        'sigma', 1.79e7);

    MatLib.titanium = struct( ...
        'eps_r', cond_eps_tensor(1.0, 2.38e6), ...
        'mu_r',  iso_tensor(1.00005), ...
        'sigma', 2.38e6);
    % [7] Magnetic materials
    MatLib.iron = struct( ...
        'eps_r', cond_eps_tensor(1.0, 1.0e7), ...
        'mu_r',  iso_tensor(4000.0), ...
        'sigma', 1.0e7);

    MatLib.nickel = struct( ...
        'eps_r', cond_eps_tensor(1.0, 1.4e7), ...
        'mu_r',  iso_tensor(600.0), ...
        'sigma', 1.4e7);

    MatLib.ferrite_nizn = struct( ...
        'eps_r', iso_tensor(15.0), ...
        'mu_r',  iso_tensor(1000.0), ...
        'sigma', 0.0);
    
    % ---------------------------------------------------------
    % 3. Preallocation
    % ---------------------------------------------------------
    num_elements = max(Mesh.ele_physTag(:, 1));

    Mesh.eps_r = repmat(MatLib.vaccum.eps_r, num_elements, 1);
    Mesh.mu_r  = repmat(MatLib.vaccum.mu_r,  num_elements, 1);

    % 추가: conductivity도 저장
    Mesh.sigma = zeros(num_elements, 1);

    Mesh.PML_cond_set = zeros(num_elements, 1);

    % ---------------------------------------------------------
    % 4. Physical Names loop
    % ---------------------------------------------------------
    physNames = Mesh.Phys.NameToID.keys();

    for i = 1:length(physNames)
        raw_name = physNames{i};
        lower_name = lower(raw_name);
        tag_id = Mesh.Phys.NameToID(raw_name);

        target_ele_ids = Mesh.ele_physTag(Mesh.ele_physTag(:, 2) == tag_id, 1);
        num_targets = length(target_ele_ids);

        if num_targets == 0
            continue;
        end

        % -----------------------------------------------------
        % 5. PML region
        % -----------------------------------------------------
        if contains(lower_name, 'pml')
            Mesh.PML_cond_set(target_ele_ids) = 1;

            % Default background material: vaccum
            bg_eps = MatLib.vaccum.eps_r;
            bg_mu  = MatLib.vaccum.mu_r;
            bg_sigma = MatLib.vaccum.sigma;

            % Find background material from physical name
            % Example: "pml_air", "pml_teflon"
            mat_keys = fieldnames(MatLib);

            for k = 1:length(mat_keys)
                key = mat_keys{k};

                if contains(lower_name, key) && ~strcmp(key, 'vaccum')
                    bg_eps = MatLib.(key).eps_r;
                    bg_mu  = MatLib.(key).mu_r;
                    bg_sigma = MatLib.(key).sigma;
                    break;
                end
            end
            num_target=length(target_ele_ids);
            if PML_set==1
                [pml_eps, pml_mu] = Assign_PML(target_ele_ids, Mesh, omega_tilde, bg_eps, bg_mu);
    
                Mesh.eps_r(target_ele_ids, :) = pml_eps;
                Mesh.mu_r(target_ele_ids, :)  = pml_mu;
                Mesh.sigma(target_ele_ids)    = bg_sigma;
                disp(['PML 적용됨']);
            elseif PML_set==2
                Mesh.eps_r(target_ele_ids, :) = repmat(bg_eps, num_target, 1);
                Mesh.mu_r(target_ele_ids, :)  = repmat(bg_mu, num_target, 1);
                Mesh.sigma(target_ele_ids)    = bg_sigma;

                disp(['[Material] PML OFF 적용됨 (배경 매질로 대체): ', physNames{i}]);

            else
                error('알 수 없는 PML_set 값입니다. (1: ON, 2: OFF)');
            end

        % -----------------------------------------------------
        % 6. Normal material region
        % -----------------------------------------------------
        else
            mat_assigned = false;
            mat_keys = fieldnames(MatLib);

            for k = 1:length(mat_keys)
                key = mat_keys{k};

                if contains(lower_name, key)
                    Mesh.eps_r(target_ele_ids, :) = ...
                        repmat(MatLib.(key).eps_r, num_targets, 1);

                    Mesh.mu_r(target_ele_ids, :) = ...
                        repmat(MatLib.(key).mu_r, num_targets, 1);

                    Mesh.sigma(target_ele_ids) = MatLib.(key).sigma;

                    mat_assigned = true;
                    break;
                end
            end

            if ~mat_assigned
                warning('Physical Tag [%s]에 해당하는 매질을 라이브러리에서 찾을 수 없어 vaccum으로 초기화합니다.', raw_name);
            end
        end
    end
end
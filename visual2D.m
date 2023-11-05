function visual2D(t_record, state_record, force_record, sep_t, param, vidFlag, vidName)
    scale = 1000;

    f = figure('position', [0, 0, 800, 800]);
    v = state_record(:,1); gamma = state_record(:,2);
    h = state_record(:,3)/scale; x = state_record(:,4)/scale;
    m = state_record(:,5); is_thrust = force_record ~= 0;

    rot = @(angle) [cos(angle) -sin(angle); sin(angle) cos(angle)];

    lv_body_len = param.l1+param.l2; lv_width = param.dia;
    lv_plf_len = param.lPLF; lv_len = param.L;
    lv_body_x = [-lv_width/2 -lv_width/2 0 lv_width/2 lv_width/2 -lv_width/2];
    lv_body_y = [-lv_len/2 -lv_len/2+lv_body_len lv_len/2 -lv_len/2+lv_body_len -lv_len/2 -lv_len/2];
    lv_body = [lv_body_y; lv_body_x]/scale;

    lv_len_2 = lv_len-param.l1; lv_width = param.dia; lv_body_len_2 = param.l2;
    lv_body_x2 = [-lv_width/2 -lv_width/2 0 lv_width/2 lv_width/2 -lv_width/2];
    lv_body_y2 = [-lv_len_2/2 -lv_len_2/2+lv_body_len_2 lv_len_2/2 -lv_len_2/2+lv_body_len_2 -lv_len_2/2 -lv_len_2/2];
    lv_body_2 = [lv_body_y2; lv_body_x2]/scale;
    
    fire_length = lv_body_len/2;
    lv_fire_x = [-lv_width/2 0 lv_width/2];
    lv_fire_y = [-lv_len/2 -lv_len/2-fire_length -lv_len/2];
    lv_fire = [lv_fire_y; lv_fire_x]/scale;

    lv_fire_x2 = [-lv_width/2 0 lv_width/2];
    lv_fire_y2 = [-lv_len_2/2 -lv_len_2/2-fire_length -lv_len_2/2];
    lv_fire_2 = [lv_fire_y2; lv_fire_x2]/scale;

    rot_lv_body = rot(pi/2-gamma(1))*lv_body;
    plot(rot_lv_body(2,:),rot_lv_body(1,:),'k','LineWidth',2);
    axis("equal"); init_limits = 5*[xlim ylim];

    ind = 1; indFFF = 1;
    while ind < length(t_record)
        lv_body_curr = lv_body;
        lv_fire_curr = lv_fire;
        if t_record(ind) > sep_t
            lv_body_curr = lv_body_2;
            lv_fire_curr = lv_fire_2;
        end

        % rotate LV body
        rot_lv_body = rot(pi/2-gamma(ind))*lv_body_curr;
        shift_lv_body = rot_lv_body+[h(ind); x(ind)].*ones(size(rot_lv_body));

        rot_lv_fire = rot(pi/2-gamma(ind))*lv_fire_curr;
        shift_lv_fire = rot_lv_fire+[h(ind); x(ind)].*ones(size(rot_lv_fire));

        if exist("h_lv_body","var") delete(h_lv_body); end
        if exist("h_lv_fire","var") delete(h_lv_fire); end
        h_lv_body = plot(shift_lv_body(2,:),shift_lv_body(1,:),'k','LineWidth',2); hold on;
        if is_thrust(ind)
            h_lv_fire = plot(shift_lv_fire(2,:),shift_lv_fire(1,:),'r','LineWidth',2);
        end

        limits = init_limits+[x(ind) x(ind) h(ind) h(ind)]; grid on;
        xlim(limits(1:2)); ylim(limits(3:4));
        title("Time since release: "+num2str(t_record(ind))+" second")
        xlabel("Downrange (km)"); ylabel("Altitude (km)");

        if vidFlag
            FFF(indFFF) = getframe(f);
            if size(FFF(indFFF).cdata,1) ~= 1800 || size(FFF(indFFF).cdata,2) ~= 1800
                ind = ind - 1; indFFF = indFFF - 1;
                continue
            end
            indFFF = indFFF + 1;
        else
            drawnow
        end
        ind = ind+1;
    end

    if vidFlag
        v = VideoWriter(vidName,"MPEG-4"); open(v);
        for i = 1:size(FFF,2)
            writeVideo(v,FFF(i));
        end
        close(v);
    end
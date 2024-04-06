function lvgraph(body)
% surface plot in body frame
upplf = body.upplf; botplf = body.botplf;
up1st = body.up1st; bot1st = body.bot1st; side1st = body.side1st;

% define coloring
cplf = [0.1 0.3 0.5]; cedgplf = [0.1 0.3 0.5];
c1st = [0.3 0.3 0.3]; cedg1st = [0 0 0];

% surface plot for the launch vehicle payload fairing
surf(upplf.x, upplf.y, upplf.z, 'FaceColor',cplf,"FaceAlpha",0.5,"EdgeColor",cedgplf); hold on;
surf(botplf.x, botplf.y, botplf.z, 'FaceColor',cplf,"FaceAlpha",0.5,"EdgeColor",cedgplf);

% surface plot for the launch vehicle body
surf(up1st.x, up1st.y, up1st.z, 'FaceColor',c1st,"FaceAlpha",0.5,"EdgeColor",cedg1st);
surf(bot1st.x, bot1st.y, bot1st.z, 'FaceColor',c1st,"FaceAlpha",0.5,"EdgeColor",cedg1st);
surf(side1st.x, side1st.y, side1st.z, 'FaceColor',c1st,"FaceAlpha",0.5,"EdgeColor",cedg1st);

xlabel("x"); ylabel("y"); zlabel("z");
axis equal; hold off

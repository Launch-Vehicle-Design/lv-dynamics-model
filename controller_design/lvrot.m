function body = lvrot(C,body,offset)
% default at rotating around the payload fairing tip
if ~exist("offset","var")
    offset = [0;0;0];
end
upplf = body.upplf; botplf = body.botplf;
up1st = body.up1st; bot1st = body.bot1st; side1st = body.side1st;

body.upplf = motion(C,upplf,offset);
body.botplf = motion(C,botplf,offset);
body.up1st = motion(C,up1st,offset);
body.bot1st = motion(C,bot1st,offset);
body.side1st = motion(C,side1st,offset);

    function update = motion(C,old,offset)
        sz = size(old.x);
        x = reshape(old.x,[1,sz(1)*sz(2)]);
        y = reshape(old.y,[1,sz(1)*sz(2)]);
        z = reshape(old.z,[1,sz(1)*sz(2)]);
        updated_xyz = C*([x; y; z]-offset);
        update.x = reshape(updated_xyz(1,:),sz);
        update.y = reshape(updated_xyz(2,:),sz);
        update.z = reshape(updated_xyz(3,:),sz);
    end

end

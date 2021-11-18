% Unit test script to test functions to convert and create line of sight
% (los) vectors within the navsu.geo namespace.


%% Test 1: calcLosEnub 
% test a simple case of [1 0 0] ecef los vector in different locations
los_xyzb = [1 0 0 1];

e_hat = [0 1 0; ... % user on ecef x axis, los vector points up
         -1 0 0];   % user on ecef y axis, los vector points west
n_hat = [0 0 1; ...
         0 0 1];
u_hat = [1 0 0; ...
         0 1 0];

% note the sign change due to flipping the direction of the vector in
% calcLosEnub
assert(all(navsu.geo.calcLosEnub(los_xyzb, e_hat, n_hat, u_hat) ...
           == [0 0 -1 1; 1 0 0 1], 'all'), ...
       'Got x-pointing los vector from various locations wrong.');

% now test multiple los vectors for user on ecef x axis
los_xyzb = [1 0 0 1; ...
            0 1 0 1; ...
            0 -1 0 1];

e_hat = [0 1 0];
n_hat = [0 0 1];
u_hat = [1 0 0];

assert(all(navsu.geo.calcLosEnub(los_xyzb, e_hat, n_hat, u_hat) ...
           == [0 0 -1 1; -1 0 0 1; 1 0 0 1], 'all'), ...
       'Got simple enub los vectors from x axis user wrong.');


% now test if opposite vectors are opposite
random_xyz = rand(5, 3);
los_xyzb = [random_xyz./vecnorm(random_xyz, 2, 2), ...
            ones(size(random_xyz, 1), 1)];
assert(all(abs(navsu.geo.calcLosEnub(los_xyzb, e_hat, n_hat, u_hat) ...
             + navsu.geo.calcLosEnub(-los_xyzb, e_hat, n_hat, u_hat)) < 1e-8, 'all'), ...
       'Opposites fail to be opposites in calcLosEnub function.');

%% Test 2: findLosEnub

% test some simple los vectors in from a position on the ecef x axis
los_xyzb = [1 0 0 1; ...
            0 1 0 1; ...
            0 -1 0 1];

e_hat = [0 1 0];
n_hat = [0 0 1];
u_hat = [1 0 0];

% note the sign change due to flipping the direction of the vector in
% calcLosEnub
assert(all(navsu.geo.findLosEnub(los_xyzb, e_hat, n_hat, u_hat) ...
           == [0 0 -1 1; -1 0 0 1; 1 0 0 1], 'all'), ...
       'Got x-pointing los vector from various locations wrong.');

% now test a large number of line of sight vectors
n_sat = 32; n_user = 10;
randLos = rand(n_sat*n_user, 3);
los_xyzb = [randLos ./ vecnorm(randLos, 2, 2), ones(n_sat*n_user, 1)];
e_hat = rand(n_user, 3);
e_hat = e_hat ./ vecnorm(e_hat, 2, 2);
n_hat = rand(n_user, 3);
n_hat = n_hat ./ vecnorm(n_hat, 2, 2);
u_hat = rand(n_user, 3);
u_hat = u_hat ./ vecnorm(u_hat, 2, 2);

% findLosEnub is essentially a call to repelem
assert(all(abs(navsu.geo.findLosEnub(los_xyzb, e_hat, n_hat, u_hat) ...
             - navsu.geo.calcLosEnub(los_xyzb, ...
                                     repelem(e_hat, n_sat, 1), ...
                                     repelem(n_hat, n_sat, 1), ...
                                     repelem(u_hat, n_sat, 1))) < 1e8, 'all'), ...
       'findLosEnub does the input rearranging wrong!');


%% Test 3: findLosXyzb

% test for two user locations, two satellite locations
xyz_user = navsu.constants.rEarth * [1 0 0; 0 1 0];
xyz_sat = (navsu.constants.rEarth + 20*1e6) * [1 0 0; 0 1 0];

los_xyzb = navsu.geo.findLosXyzb(xyz_user, xyz_sat);

% the first and last should be directly overhead
assert(all(los_xyzb([1 4], 1:3) == [1 0 0; 0 1 0], 'all'), ...
       'Satellites are not overhead where they should be!');

assert(all(diag(los_xyzb([2 3], [1 2])) < 0, 'all'), ...
       'Satellite offset from vertical in wrong direction!');

assert(los_xyzb(2, 1) == los_xyzb(3, 2) && los_xyzb(3, 1) == los_xyzb(2, 2), ...
       'Symmetric los vectors are not symmetric :-(');



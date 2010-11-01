#### Congfiguration ####

require 'catpaws/ec2'
require 'pathname'

set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :key, ENV['EC2_KEY'] #this should be the name of the key in aws, not the actual keyfile
set :ssh_options, {
  :keys => [ENV['EC2_KEYFILE']],
  :user => "ubuntu"
}
set :ami, 'ami-52794c26' #32-bit ubuntu lucid server (eu-west-1)
set :instance_type, 'm1.small'
set :working_dir, '/mnt/work'

set :group_name, 'role_of_rest_in_astrocytes'
set :nhosts, 1

set :snap_id, `cat SNAPID`.chomp #empty until you've created a snapshot
set :vol_id, `cat VOLUMEID`.chomp #empty until you've created a new volume
set :ebs_size, 5  
set :availability_zone, 'eu-west-1a'  #wherever your ami is. 
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'

set :git_url, 'http://github.com/cassj/role_of_rest_in_astrocytes/raw/master'

#define the mount points of subprojects. Note that these are also defined in the 
#subproject Capfiles, so you need to make sure they agree if you change them.
set :sp_johnson, '/mnt/johnson'
set :sp_rest_chip_astro, '/mnt/rest_chip_astro'
set :sp_dnrest_astro_xpn, '/mnt/dnrest_astro_xpn'


# Try and load a local config file to override any of the above values, should one exist.
# So that if you change these values, they don't get overwritten if you update the repos.
begin
 load("Capfile.local")
rescue Exception
end

#cap EC2:start
#cap EBS:create (unless you want to use the one that already exists)
#cap EBS:attach
#cap EBS:format_xfs (unless you're using a snapshot)
#cap EBS:mount_xfs
#
# For each subproject that you want to use, 
# cd subproject_dir
# cap EBS:create
# cap EBS:attach
# cap EBS:mount_xfs
# cd ../
#
# run tasks as required
# 
#cap EBS:snapshot
#cap EBS:unmount
#cap EBS:detach
#cap EBS:delete
#cap EC2:stop



### Instance Setup Tasks
desc "install R on all running instances in group group_name"
task :install_r, :roles  => group_name do
  user = variables[:ssh_options][:user]
  run "mkdir -p #{working_dir}/scripts"
  sudo 'apt-get update'
  sudo 'apt-get -y install r-base'
  sudo 'apt-get -y install build-essential libxml2 libxml2-dev libcurl3 libcurl4-openssl-dev'
  run "cd #{working_dir}/scripts &&  curl '#{git_url}/scripts/R_setup.R' > R_setup.R"
  run "chmod +x #{working_dir}/scripts/R_setup.R"
  sudo "Rscript #{working_dir}/script/R_setup.R"
end
before "install_r", "EC2:start"
  

#### Tasks ####
desc "overlap ESC NSC and Astro sites"
task :esc_nsc_astro_overlap, :roles => group_name do
   run "mkdir -p #{working_dir}/scripts"
   run "cd #{working_dir}/scripts && curl #{git_url}/overlap_nsc_esc_astro.rnw overlap_nsc_esc_astro.rnw"  
end
before "esc_nsc_astro_overlap", "EC2:start"


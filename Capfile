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


#should be using the buckley_lab_ami but it isn't working yet.
set :ami, 'ami-7e5c690a' #32-bit Ubuntu Maverick, eu-west-1
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



### Instance Setup Tasks.

#these will be done on buckley ami when I've got it working  

desc "install latest R-devel"
task :install_r, :roles  => group_name do
  user = variables[:ssh_options][:user]
#  sudo "apt-get -y build-dep r-base"  #everything for R
#  sudo "apt-get -y install libxml2 libxml2-dev libcurl3 libcurl4-openssl-dev" #and for some of the Bioconductor stuff
	
  run "mkdir -p #{working_dir}/scripts"
  upload("scripts/R_setup.R","#{working_dir}/scripts/R_setup.R")
  run "cd #{working_dir}/scripts && chmod +x R_setup.R"

  r_src_url = 'ftp://ftp.stat.math.ethz.ch/Software/R/R-devel.tar.bz2'
  run "cd #{working_dir} && curl #{r_src_url} > R-devel.tar.bz2"
  run "cd #{working_dir} && tar -xvjf  R-devel.tar.bz2"
  run "cd #{working_dir}/R-devel && umask 022 && ./configure"  
  run "cd #{working_dir}/R-devel && umask 022 && make"
  run "cd #{working_dir}/R-devel && umask 022 && sudo make install"

  sudo "Rscript #{working_dir}/scripts/R_setup.R"

end
before "install_r", "EC2:start"
  



#### Tasks ####

#add files in the form 'remote_src_path' => 'local_target'
dl_files = Hash.new()

desc "overlap ESC NSC and Astro sites"
task :overlap_nsc_esc_astro, :roles => group_name do
  run "mkdir -p #{working_dir}/scripts"
  upload('scripts/overlap_nsc_esc_astro.R', "#{working_dir}/scripts/overlap_nsc_esc_astro.R")
  run "chmod +x #{working_dir}/scripts/overlap_nsc_esc_astro.R"
  run "cd #{mount_point} && Rscript #{working_dir}/scripts/overlap_nsc_esc_astro.R #{sp_johnson}/NS5/PET/RangedData.R #{sp_johnson}/ESC/PET/RangedData.R #{sp_rest_chip_astro}/Macs/NA_peaks.RangedData.RData #{mount_point}/overlap_nsc_esc_astro"
end
before "overlap_nsc_esc_astro", "EC2:start"


desc "fetch_results"
task :fetch_results, :roles => group_name do
  files = capture("ls #{mount_point}/nsc_esc_astro_overlap_*").split("\n")
  files.each {|f|
    fname = f.sub(/^.*\//,'')
    download(f,"results/#{fname}" )
  }
end 
before "fetch_results", "EC2:start"



# tests_install

# 1/ Install GENBO packages

# 2/ Launch Tests Paths
perl -I $GENBO/lib/obj-nodb/ tests_paths.pl

# 3/ Launch Tests PolyQuery Interface:
perl -I $GENBO/lib/obj-nodb/ tests_polyquery.pl -path_polyquery=$polymorphism-cgi/json_output_nodb/


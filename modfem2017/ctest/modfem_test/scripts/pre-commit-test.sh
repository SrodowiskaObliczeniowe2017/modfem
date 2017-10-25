# Run fast test to check
## if all selected targets are building properly
## if fast funcional tests are working ok
mf_root_path="$(git rev-parse --show-toplevel)"
mf_test_dirs="$(ls -d $mf_root_path/bin_cmake/*test*/)"

for mf_test_dir in $mf_test_dirs; do
	echo "Testing $mf_test_dir"

	if [ -d $mf_test_dir ]; then
		echo "$(cd $mf_test_dir)"
	else
		echo "Cannot access directory $mf_test_dir"
		exit 1
	fi

	if [ ! "$(ctest -R Make)" ]; then
		echo "Building OK."
	else
		echo "Building failed."
		exit 1
	fi

	if [ ! "$(ctest -R test)" ]; then
		echo "Functional tests OK."
	else
		echo "Functional tests failed."
		exit 1
	fi
done

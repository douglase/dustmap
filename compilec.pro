pro compilec, routine_name, routine_dir

if n_elements(routine_dir) eq 0 then routine_dir = './'

make_dll, routine_name, routine_name, routine_name, input=routine_dir, output=routine_dir,cc='gcc -c -fPIC -O3 %C -o %O'

end

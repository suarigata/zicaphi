#source ~/xeon-phi-offload-example/configure.sh
source ~edson/src/cwp/set_path.sh
diff c.su ../cmp/c.su
diff cmp.stack.su ../cmp/cmp.stack.su
diff cmp.coher.su ../cmp/cmp.coher.su
sucmp cmp.stack.su ../cmp/cmp.stack.su && sucmp cmp.coher.su ../cmp/cmp.coher.su && sucmp c.su ../cmp/c.su 


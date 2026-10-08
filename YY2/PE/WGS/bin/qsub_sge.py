# Patched by Codex for YY2 WGS: fixed retry resources and bounded qw wait.
#!/mnt/gpfs/Users/wangning/software/Miniconda/bin/python3
########################################## import ################################################
import argparse, os, sys, logging, glob, sqlite3, socket, subprocess, io, time, re
from datetime import datetime
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append('/mnt/gpfs/Users/wangning/program/sge')
import mysqlite
############################################ ___ #################################################
__doc__ = ''
__author__ = ''
__mail__ = ''
__date__ = ''
__version__ = '1.0.0'
############################################ main ##################################################
def mylogger(log_file):
	logger = logging.getLogger('mylogger')
	logger.setLevel(logging.DEBUG)
	fh = logging.FileHandler(log_file)
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	return(logger)

def report(level,info):
	date_now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	if level == "ERROR":
		sys.stderr.write("{0} - {1} - ERROR - {2}\n".format(date_now,os.path.basename(__file__),info))
		sys.exit(1)
	elif level == "INFO":
		sys.stdout.write("{0} - {1} - INFO - {2}\n".format(date_now,os.path.basename(__file__),info))
	elif level == "DEBUG":
		sys.stdout.write("{0} - {1} - DEBUG - {2}\n".format(date_now,os.path.basename(__file__),info))
		sys.exit(1)
	return()

def check_file(file):
	if os.path.exists(file):
		return(os.path.abspath(file))
	else:
		info = "{0} does not exis!".format(file)
		report("ERROR",info)

def check_dir(dir):
	dir = os.path.abspath(dir)
	if not os.path.exists(dir):
		os.system("mkdir -p {0}".format(dir))
		log.info("mkdir {0}".format(dir))
	return(dir)

def run_cmd(cmd):
	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=-1)
	proc.wait()
	stream_stdout = io.TextIOWrapper(proc.stdout, encoding='utf-8')
	stream_stderr = io.TextIOWrapper(proc.stderr, encoding='utf-8')
	stdout = str(stream_stdout.read())
	stderr = str(stream_stderr.read())
	return(stdout, stderr)

def create_database(db):
	conn = sqlite3.connect(db)
	cmd = '''CREATE TABLE JOB
		(SCRIPT   TEXT    PRIMARY KEY ,
		 QUEUE    TEXT ,
		 P        TEXT ,
		 VF       TEXT ,
		 HVMEM    TEXT ,
		 MAXCYCLE TEXT ,
		 JOBID    TEXT ,
		 TRYTIMES TEXT ,
		 STATE    TEXT) ; '''
	conn.execute(cmd)
	conn.commit()

def split_file(shell, line_num, prefix, qsub_dir):
	n = 1
	line_dic = {}
	sub_script_list = []
	with open(shell, 'r') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			if n in line_dic:
				if len(line_dic[n]) == line_num:
					n += 1
				line_dic.setdefault(n, []).append(line)
			else:
				line_dic.setdefault(n, []).append(line)
	for m, comm_list in line_dic.items():
		sub_script = '{}/{}_{}.sh'.format(qsub_dir, prefix, str(m).zfill(4))
		sub_script_list.append(sub_script)
		comm_list.append('echo This-Work-is-Completed!\n')
		with open(sub_script, 'w') as OUT:
			OUT.write(' && \\\n'.join(comm_list))
	return(sub_script_list)

def qsub_job(run_jobs, wait_jobs):
	for idx in range(maxjob - len(run_jobs)):
		try:
			sub_script, queue, pe, vf, hvmem, maxcycle, jobid, trytimes, state = list(wait_jobs[idx])
			maxcycle = int(maxcycle)
			trytimes = int(trytimes)
			# Keep resource requests fixed across retries.  The historical
			# exponential bump can strand jobs in qw when the cluster has
			# enough ordinary slots but no node matching the inflated memory.
			vf = round(float(vf), 2)
			hvmem = round(float(hvmem), 2)
			sub_script_dir = os.path.dirname(sub_script)
			cmd = 'cd {sub_script_dir} && qsub -cwd -pe smp {pe} -l vf={vf}G,h_vmem={hvmem}G -q {queue} {sub_script}'.format(sub_script_dir = sub_script_dir, pe = pe, vf = vf, hvmem = hvmem, queue = queue, sub_script = sub_script)
			stdout, stderr = run_cmd(cmd)
			trytimes = int(trytimes) + 1
			log.info('{sub_script} try {trytimes} times with parameters: -pe smp {pe} -l vf={vf}G,h_vmem={hvmem}G -q {queue}'.format(sub_script = sub_script, trytimes= trytimes, pe = pe, vf = vf, hvmem = hvmem, queue = queue))
			if stderr:
				jobid = '-'
				state = 'waiting'
				log.error('qsub {sub_script} failed with error: {stderr}'.format(sub_script = sub_script, stderr = stderr))
			else:
				parts = stdout.split()
				if parts and parts[0].isdigit():
					jobid = parts[0]
				else:
					jobid = parts[2]
				state = 'running'
				log.info(stdout.split('\n')[0])
			db.update({'SCRIPT':sub_script}, {'VF':vf, 'HVMEM': hvmem, 'JOBID':jobid, 'TRYTIMES':trytimes, 'STATE':state})
		except IndexError as e:pass
		
def qacct_job(run_jobs):
	for run_job in run_jobs:
		sub_script, queue, pe, vf, hvmem, maxcycle, jobid, trytimes, state = list(run_job)
		stdout, stderr = run_cmd('qacct -j {}'.format(jobid))
		if stderr:
			stdout, stderr = run_cmd('qstat -j {}'.format(jobid))
			if not stderr and max_queue_minutes > 0:
				state_out, state_err = run_cmd("qstat -u $USER -s a | awk '$1==\"{}\"{{print $5}}'".format(jobid))
				if state_out.strip() == 'qw':
					submitted = None
					for line in stdout.split('\n'):
						if line.startswith('submission_time:'):
							value = line.split(':', 1)[1].strip()
							try:
								submitted = datetime.strptime(value, '%a %b %d %H:%M:%S %Y')
							except Exception:
								submitted = None
					if submitted and (datetime.now() - submitted).total_seconds() > max_queue_minutes * 60:
						qdel_out, qdel_err = run_cmd('qdel {}'.format(jobid))
						log.warning('{} stayed qw longer than {} minutes; qdel and reset to waiting without resource escalation'.format(jobid, max_queue_minutes))
						db.update({'SCRIPT':sub_script}, {'JOBID':'-', 'TRYTIMES':'0', 'STATE':'waiting'})
						continue
			if stderr.startswith('Following jobs do not exist'):
				time.sleep(30)
				stdout, stderr = run_cmd('qacct -j {}'.format(jobid))
				if stderr:
					state = 'waiting'
				else:
					for line in stdout.split('\n'):
						if line.startswith('exit_status'):
							line = line.split()
							if line[1] == '0':
								state = 'done'
							else:
								if int(trytimes) == int(maxcycle):
									state = 'error'
								else:
									state = 'waiting'
				db.update({'SCRIPT':sub_script}, {'STATE':state})
				print('数据库更新为waiting，该处查询可能有问题，请注意')
			continue
		for line in stdout.split('\n'):
			if line.startswith('exit_status'):
				line = line.split()
				if line[1] == '0':
					state = 'done'
				else:
					if int(trytimes) == int(maxcycle):
						state = 'error'
					else:
						state = 'waiting'
				db.update({'SCRIPT':sub_script}, {'STATE':state})

def jobs_guide():
	run_jobs = db.ask({'STATE':'running'})
	wait_jobs = db.ask({'STATE':'waiting'})
	error_jobs = db.ask({'STATE':'error'})
	done_jobs = db.ask({'STATE':'done'})
	# 查询运行中的任务状态
	if run_jobs:
		qacct_job(run_jobs)
	if wait_jobs:
		qsub_job(run_jobs, wait_jobs)
	if not wait_jobs and not run_jobs:
		if error_jobs:
			for i in error_jobs:
				log.error('{} run failed!'.format(i[0]))
#				report('INFO','{} run failed!'.format(i[0]))
			report('ERROR', '{} run failed!'.format(i))
		else:
			log.info('all jobs run successfully!')
		return(True)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument( help='input file', dest='input')
	parser.add_argument('-m', '--mem', help='request memory (G), default is [1]', dest='mem', type=int, default=1)
	parser.add_argument('-t', '--thread', help='request slot range for parallel jobs, default is [1]', dest='pe', type=int, default=1)
	parser.add_argument('-q', '--queue', help='job queue, default is [all.q]', dest='queue', type=str, default='all.q')
	parser.add_argument('-l', '--line', help='line number for sub_script, default is [1]', dest='line', type=int, default=1)
	parser.add_argument('-j', '--maxjob', help='max job number for qsub, default is [30]', dest='maxjob', type=int, default=30)
	parser.add_argument('-c', '--maxcycle', help='number of times for qsub, default is [1]', dest='cycle', type=int, default=1)
	parser.add_argument('-p', '--prefix', help='prefix for job, do not start with digits, default is script name', dest='prefix', type=str, required=False)
	parser.add_argument('--max-queue-minutes', help='cancel and requeue jobs that remain qw longer than this many minutes; 0 disables this guard, default [360]', dest='max_queue_minutes', type=int, default=360)
	args=parser.parse_args()
	report("INFO", "Start")
	# 检查当前主机
	if socket.gethostname() != 'mgt':
		report('ERROR', 'The current host is \"{}\", please ssh to the \"master\" node and try again!'.format(socket.gethostname()))
	# 定义文件
	global log, maxjob
	shell = check_file(args.input)
	shell_db = shell + '.db'
	log = mylogger(shell + '.log')
	maxjob = args.maxjob
	global max_queue_minutes
	max_queue_minutes = args.max_queue_minutes
	qsub_dir = check_dir(shell + '.qsub')
	if args.prefix:
		prefix = args.prefix
	else:
		prefix = os.path.splitext(shell.split('/')[-1])[0].lstrip('0123456789_-.')
	# 防止以数字开头
	if re.compile(r"^[0-9].*").match(prefix):
		prefix = 'oe_' + prefix
	# 检查创建数据库
	if not os.path.exists(shell_db):
		create_database(shell_db)
		log.info('create database: {}'.format(shell_db))
		# 拆分脚本
		sub_script_list = split_file(shell, args.line, prefix, qsub_dir)
		log.info('split {shell} to {num} subscripts with prefix: {prefix}'.format(shell = shell, num = len(sub_script_list), prefix = prefix))
		conn = sqlite3.connect(shell_db)
		for sub_script in sub_script_list:
			cmd = 'INSERT INTO JOB VALUES(\"{SCRIPT}\", \"{QUEUE}\", \"{P}\", \"{VF}\", \"{HVMEM}\", \"{MAXCYCLE}\", \"{JOBID}\", \"{TRYTIMES}\", \"{STATE}\");'.format(SCRIPT = sub_script, QUEUE = args.queue, P = args.pe, VF = args.mem, HVMEM = args.mem, MAXCYCLE = args.cycle, JOBID = "-", TRYTIMES = "0", STATE = "waiting")
			conn.execute(cmd)
		conn.commit()
	global db
	db = mysqlite.mydb(shell_db, 'JOB')
	# 更新数据库中失败任务信息，准备重新投递
	error_jobs = db.ask({'STATE':'error'})
	if error_jobs:
		for error_job in error_jobs:
			sub_script = list(error_job)[0]
			db.update({'SCRIPT':sub_script}, {'VF':args.mem, 'HVMEM': args.mem, 'JOBID':'-', 'TRYTIMES':'0', 'STATE':"waiting"})
	run_jobs = db.ask({'STATE':'running'})
	# 更新数据库中正在运行中的任务信息，准备重新投递
	if run_jobs:
		qacct_job(run_jobs)
	# 循环监控
	while True:
		s = jobs_guide()
		if s == None:pass
		else:break
		log.info('sleep 120 seconds')
		time.sleep(120)
		log.info('wake up')

if __name__=="__main__":
	main()
	report("INFO", "End")

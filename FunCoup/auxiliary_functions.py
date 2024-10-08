import joblib
import contextlib
import zipfile
import gzip
import os 
import requests, json

def downloadFile(urlToDownload, filename):
    # Downloading file from url
    url = urlToDownload
    myfile = requests.get(url)
    if myfile.status_code==200:
        open(filename, "wb").write(myfile.content)
    else:
        print("Cant download file from: "+url)

def unzipFile(filename,path):
    # Unzipping file to path
    print("unzipping file: "+filename+".zip")
    with zipfile.ZipFile(path+filename+".zip", "r") as zip_ref:
       zip_ref.extractall( path=path)

def gzipFile(filename,path):
    # Unzipping file to path
    print("unzipping file: "+filename+".gz")
    with gzip.open(path+filename+".gz", "r") as gzip_ref:
       gzip_content = gzip_ref.read()
    with open(path+filename, 'wb') as fout:
        fout.write(gzip_content)

def downloadAndUnzipFile(url,format='zip'):
    if not os.path.exists("data/tmp/"+url.split("/")[-1].split(".zip")[0]): # Might already have been downloaded.
        downloadFile(url, "data/tmp/"+url.split("/")[-1])
    if format == 'zip':
        if os.path.exists("data/tmp/"+url.split("/")[-1]):
            unzipFile(url.split("/")[-1].split(".zip")[0],"data/tmp/")
    elif format == 'gz':
        if os.path.exists("data/tmp/"+url.split("/")[-1]):
            gzipFile(url.split("/")[-1].split(".gz")[0],"data/tmp/")

def removeTempFile(filename):
    # Removing file named filename
    if os.path.isfile(filename):
        os.remove(filename)

def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

def calculateChunkIndex(len_df,num_cores):
    # Calculate the chunk size
    chunk_size = len_df // num_cores
    # Initialize a list to store start and end indices for each chunk
    chunk_indices = []
    # Create chunks and calculate start and end indices
    for i in range(num_cores):
        start_index = i * chunk_size
        end_index = (i + 1) * chunk_size if i < num_cores - 1 else len_df
        chunk_indices.append((start_index, end_index))
    
    return chunk_indices

def sortPair(a,b):
    if a>b: return a,b
    else: return b,a

def sortPair2(a,b):
    tmp_a = int(a.split('_')[0])
    tmp_b = int(b.split('_')[0])
    if tmp_a>tmp_b: return a,b
    else: return b,a

def print_config(config, char='#', length=40):
    frame = char * length
    print(frame)
    for key,val in config['instance'].items():
        print('%s: %s' %(key,','.join(val)))
    print(frame)

@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()

def print_frame(message, char='#', length=40):
    """Prints a frame with a message in the center.

    Args:
        message (str): The message to display in the frame.
        char (str, optional): The character used to create the frame. Default is '#'.
        length (int, optional): The length of the frame. Default is 40.
    """
    frame = char * length
    centered_message = message.center(length - 2)
    print(f"{frame}\n{char} {centered_message} {char}\n{frame}")

def create_log_file(file_path,content):
    with open(file_path, "w") as f:
        f.write(content)
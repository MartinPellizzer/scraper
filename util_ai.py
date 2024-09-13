from llama_cpp import Llama

llm = None
llm_loaded = False

vault = '/home/ubuntu/vault'
llms_folderpath = f'{vault}/llms'

model_path='/home/ubuntu/vault/llms/Meta-Llama-3-8B-Instruct.Q8_0.gguf'
model_path='/home/ubuntu/vault/llms/Meta-Llama-3.1-8B-Instruct-Q4_K_M.gguf'

def gen_reply(prompt, model=None):
    global llm
    global llm_loaded
    global model_path
    if model:
        model_path = f'{llms_folderpath}/{model}'
    print('************************************')
    print(model_path)
    print('************************************')
    n_ctx = 8192
    n_ctx = 20000
    if llm_loaded == False:
        llm = Llama(
            model_path=model_path,
            n_gpu_layers=-1,
            n_ctx = n_ctx,
        )
        llm_loaded = True

    stream = llm.create_chat_completion(
        messages = [{
            'role': 'user',
            'content': prompt,
        }],
        stream = True,
        temperature = 1.0,
    )
    
    reply = ''
    for chunk in stream:
        if 'content' in chunk['choices'][0]['delta'].keys():
            token = chunk['choices'][0]['delta']['content']
            reply += token
            print(token, end='', flush=True)

    return reply


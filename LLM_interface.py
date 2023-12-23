from PIL import Image
import pandas as pd
import requests

from openai_key import *
from openai import OpenAI

def textgen_api(prompt, url = 'http://localhost:5000/api/v1/generate'):
  request = {
      'prompt': prompt,
      'max_new_tokens': 500,
      'auto_max_new_tokens': False,
      'max_tokens_second': 0,

      # Generation params. If 'preset' is set to different than 'None', the values
      # in presets/preset-name.yaml are used instead of the individual numbers.
      'preset': 'None',
      'do_sample': True,
      'temperature': 0.7,
      'top_p': 0.1,
      'typical_p': 1,
      'epsilon_cutoff': 0,  # In units of 1e-4
      'eta_cutoff': 0,  # In units of 1e-4
      'tfs': 1,
      'top_a': 0,
      'repetition_penalty': 1.18,
      'repetition_penalty_range': 0,
      'top_k': 40,
      'min_length': 0,
      'no_repeat_ngram_size': 0,
      'num_beams': 1,
      'penalty_alpha': 0,
      'length_penalty': 1,
      'early_stopping': False,
      'mirostat_mode': 0,
      'mirostat_tau': 5,
      'mirostat_eta': 0.1,
      'grammar_string': '',
      'guidance_scale': 1,
      'negative_prompt': '',

      'seed': -1,
      'add_bos_token': True,
      'truncation_length': 4000,
      'ban_eos_token': False,
      'custom_token_bans': '',
      'skip_special_tokens': True,
      'stopping_strings': []
  }
  response = requests.post(url, json=request)
#   print(response.status_code)
  if response.status_code == 200:
      result = response.json()['results'][0]['text']
  return result

def gpt_3_turbo_completions(prompt):
   client = OpenAI(api_key=gpt_key,)
   response = client.completions.create(
    model = "gpt-3.5-turbo-instruct",
    prompt = prompt,
    max_tokens = 150
    )

   return response.choices[0].text.strip()


def gpt_3_turbo_chat(prompt):
    client = OpenAI(api_key=gpt_key,)
    response = client.chat.completions.create(
        messages=[
            {
                "role": "user",
                "content": prompt,
            }
        ],
        model = "gpt-3.5-turbo",
        max_tokens = 150
        )

    return response.choices[0].message.content

def gpt_4_turbo_chat(prompt):
    client = OpenAI(api_key=gpt_key,)
    response = client.chat.completions.create(
        messages=[
            {
                "role": "user",
                "content": prompt,
            }
        ],
        model = "gpt-4",
        max_tokens = 150
        )
    # print(response.choices[0])
    # print(response.choices[0].message.content)
    return response.choices[0].message.content
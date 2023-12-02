from LLM_interface import *
def analysis_introduction(introduction):
    prompt = introduction +'\n\n'
    prompt += 'Task: The text above is an introduction to a text. '
    prompt += 'The text above is an introduction to a text. I now want to find the reference mentioned in it about the growth of certain cells in different oxygen concentration environments.\n'
    prompt += 'I hope you will first analyze whether each reference is related to what I need, and then give a Yes/No conclusion to each reference.\n'

    prompt += 'Output Format Example:\n'
    prompt += '[<citation number>] <analyze> | Coclusion:Yes/No.\n'
    # prompt += 'Please note that <reference index number in the paper> refers to the citation number that appears in the article.'


    result = gpt_3_turbo_chat(prompt)
    return result

def analysis_citation(citation):
    prompt = citation +'\n\n'
    prompt += 'Task: The above paragraph is quoted from an article.'
    prompt += 'Is this article about the growth of certain cells in different oxygen concentration environments?\n'
    prompt += 'Please analyze first before giving a yes/no conclusion.\n'

    prompt += 'Output Format Example:\n'
    prompt += ' Analysis process | Coclusion:Yes/No.\n'
    # prompt += 'Please note that <reference index number in the paper> refers to the citation number that appears in the article.'


    result = gpt_3_turbo_chat(prompt)
    return result
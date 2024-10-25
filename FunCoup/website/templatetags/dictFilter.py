from django import template

register = template.Library()

# Here are some templatetags, used to simplify collection of information from the django-html integration.
# All of these are used in the html of the website

# To get value from any dict using key
@register.simple_tag
def get_value(dictionary, key):
    return dictionary.get(key)

# To get link info from the linkDict, where keys are a combination of nodeIDa|nodeIDb
@register.simple_tag
def get_link(dictionary, nodeA, nodeB):
    if nodeA+"|"+nodeB in dictionary:
        return dictionary[nodeA+"|"+nodeB]
    else:
        return dictionary[nodeB+"|"+nodeA]
    
# To get information from a list at a given index
@register.filter
def index(sequence, position):
    return sequence[position]

# To get a pathway name from the pathway dict using a pathway ID
@register.simple_tag
def get_pathway(dictionary, key):
    return dictionary.get(key).split("(pathwayID")[0].strip()

import logging
#logging.basicConfig(level = logging.DEBUG)
## Saving it to a logging file
logging.basicConfig(filename = 'test.log', level = logging.DEBUG, format = '%(asctime)s:%(levelname)s:')
## Changing the format




def add(x,y):
    return x + y

num_1 = 10
num_2 = 5 
add_result = add(num_1, num_2)
logging.debug('Add: {} + {} = {}'.format(num_2, num_2, add_result))

# 3rd party packages
# parseConfigs=['pkg-config --cflags --libs playerc++']

parseConfigs=['pkg-config --cflags --libs playerc++']
env = Environment (
  CC = 'g++',
  CCFLAGS = Split ('-pthread -pipe  -W -Wall -O2'),
)


# Parse all the pacakge configurations
for cfg in parseConfigs:
  env.ParseConfig(cfg)

env.Program('main','main.cpp',)

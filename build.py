import os
import shutil

from photometrus.settings import gen_user_prime_dir, gen_config_dir, gen_pipeline_file_name

# sex.conv -> abs/path/sex.conv
# tempsource.param -> abs/path/tempsource.param
# default.nnw


def copy_configs(overwrite=False):
    if not os.path.isdir(gen_config_dir()) or overwrite:
        shutil.copytree(os.path.join(gen_pipeline_file_name(), 'configs'), gen_config_dir(), dirs_exist_ok=True)


def update_configs(config_dir, abs_path_list=('sex.conv', 'tempsource.param', 'default.nnw')):
    configs = [os.path.join(config_dir, f) for f in os.listdir(config_dir) if f.endswith('.config')]
    for config in configs:
        with open(config, 'r') as f:
            config_txt = f.read()
        for abs_path in abs_path_list:
            config_txt = config_txt.replace(abs_path, os.path.join(config_dir, abs_path))
        with open(config, 'w') as f:
            f.write(config_txt)


def main():
    config_dir = gen_config_dir()
    copy_configs(overwrite=True)
    update_configs(config_dir)


if __name__ == '__main__':
    main()

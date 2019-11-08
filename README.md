# sp-pharkin

Library that can be used for basic PharmacoKinetics calculations.

I wrote this code while studying basic pharmacokinetiks theory from the book:
"Pharmacokinetics" by Philip Rowe, and incoprorates formulas and practice questions in the form of tests.
## Main Features
Work in progress!!
<!---
- Simulation enviroment follows [OpenAI gym](https://github.com/openai/gym) and [rllab](https://github.com/rll/rllab) APIs. It returns observation, reward, done, info at each step, which means the simulator is "reinforcement-learning-ready".
-->



## Installation
<!---
[comment]: It is highly recommended to use `pip` to install `simglucose`, follow this [link](https://pip.pypa.io/en/stable/installing/) to install pip.

Auto installation:
```bash
pip install simglucose
```
-->


Manual installation: 
```bash
git clone https://github.com/spitoglou/sp-pharkin.git
cd sp-pharkin
```
If you have `pip` installed, then
```bash
pip install -e .
```
If you do not have `pip`, then
```bash
python setup.py install
```


## Reporting issues
Report any bugs, enhancements or even discussion by [creating issues]https://github.com/spitoglou/sp-pharkin/issues/new).

## How to contribute
The following instruction is originally from the [contribution instructions of sklearn](https://github.com/scikit-learn/scikit-learn/blob/master/CONTRIBUTING.md).

The preferred workflow for contributing to simglucose is to fork the
[main repository](https://github.com/spitoglou/sp-pharkin) on
GitHub, clone, and develop on a branch. Steps:

1. Fork the [project repository](https://github.com/spitoglou/sp-pharkin)
   by clicking on the 'Fork' button near the top right of the page. This creates
   a copy of the code under your GitHub user account. For more details on
   how to fork a repository see [this guide](https://help.github.com/articles/fork-a-repo/).

2. Clone your fork of the simglucose repo from your GitHub account to your local disk:

   ```bash
   $ git clone git@github.com:YourLogin/sp-pharkin.git
   $ cd sp-pharkin
   ```

3. Create a ``feature`` branch to hold your development changes:

   ```bash
   $ git checkout -b my-feature
   ```

   Always use a ``feature`` branch. It's good practice to **never work on the ``master`` branch**!

4. Develop the feature on your feature branch. Add changed files using ``git add`` and then ``git commit`` files:

   ```bash
   $ git add modified_files
   $ git commit
   ```

   to record your changes in Git, then push the changes to your GitHub account with:

   ```bash
   $ git push -u origin my-feature
   ```

5. Follow [these instructions](https://help.github.com/articles/creating-a-pull-request-from-a-fork)
to create a pull request from your fork. This will send an email to the committers.

(If any of the above seems like magic to you, please look up the
[Git documentation](https://git-scm.com/documentation) on the web, or ask a friend or another contributor for help.)

## License

sp-pharkin is licensed under the [MIT License](http://opensource.org/licenses/MIT).

Copyright 2015 **Stavros Pitoglou**
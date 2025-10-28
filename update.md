Make new release
================

- Update revision in `eos/__init__.py`
- Commit changes `git commit -a -m "your message here"`
- Tag version `git tag v3.x.y`
- Push changes `git push` and `git push --tags`
- This should trigger the **Release** action on GitHub that builds a new version and uploads it to PyPI.


Update on AMOR
==============

- Login via SSH using the **amor** user.
- Activate eos virtual environment `source /home/software/virtualenv/eosenv/bin/activate`
- Update eos packge `pip install --upgrade amor-eos`

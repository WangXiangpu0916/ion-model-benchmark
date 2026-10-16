# 项目工作约定

## 版本控制

- 本项目使用私密仓库 https://github.com/WangXiangpu0916/Ion-Model-Benchmark，远程为 `origin`，主分支为 `main`。
- 用户已授权：以后每次完成本项目的修改，运行与改动相关的检查；全部通过后，自动提交本次修改并推送到此私密仓库，无需再次询问。
- 提交前检查 diff，仅提交本次任务涉及的修改；遵守现有 `.gitignore`，保留仓库的私密状态。
- 检查失败或发现冲突时，先修复并复查；无法解决时说明阻碍，不将未通过检查的修改推送。
- 推送后确认本地提交与远程分支同步，并在回复中简要报告检查结果和提交编号。

## Windows shell

- 优先使用 `pwsh.exe -NoProfile`；只有需要 Windows PowerShell 5.1 兼容性时才使用 `powershell.exe`。

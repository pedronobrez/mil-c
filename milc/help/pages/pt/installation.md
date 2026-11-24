---
title: Instalação
---

O MIL-C é publicado como instaladores na página de releases do repositório,
`https://github.com/pedronobrez/mil-c/releases`. Cada versão tem cinco arquivos, um ou dois por
sistema, e um arquivo `SHA256SUMS` para conferi-los. Cada instalador traz o seu próprio Python,
Qt, numpy e código numérico compilado: não é preciso instalar mais nada antes, e nenhum deles
precisa de R.

## Os cinco instaladores

Para a versão 1.0.0:

| Arquivo | Sistema | O que faz |
|---|---|---|
| `MIL-C-1.0.0-macos-arm64.dmg` | macOS em Apple Silicon | uma imagem de disco: arraste **MIL-C** para **Applications** |
| `MIL-C-1.0.0-windows-x64.msi` | Windows, 64 bits | instala para o usuário atual, com uma entrada no menu Iniciar |
| `MIL-C-1.0.0-windows-x64.zip` | Windows, 64 bits | o mesmo programa como pasta: descompacte onde quiser e rode `MIL-C.exe` |
| `MIL-C-1.0.0-linux-x86_64.deb` | Debian, Ubuntu, Mint | instala em `/opt/MIL-C`, com entrada no menu e `milc` no PATH |
| `MIL-C-1.0.0-linux-x86_64.tar.gz` | qualquer Linux, x86-64 | descompacte e rode `install.sh`, ou rode o programa onde está |

Observações sobre cada um:

- **macOS.** Só Apple Silicon: o numba, com que o UMAP e o FlowSOM são compilados, não publica
  mais versões para Macs Intel.
- **Windows `.msi`.** Instala na sua própria pasta `AppData\Local\MIL-C` e não pede
  administrador. Uma versão mais nova substitui a anterior em vez de ficar ao lado dela. Para
  remover, use **Configurações ▸ Aplicativos** como para qualquer programa.
- **Windows `.zip`.** Para uma máquina em que não se pode rodar instalador. Nada é gravado fora da
  pasta onde você o descompacta, a não ser as preferências e os caches do programa.
- **Linux `.deb`.** `sudo apt install ./MIL-C-1.0.0-linux-x86_64.deb`. O pacote declara as
  poucas bibliotecas do sistema de que o Qt precisa, e o apt as instala junto. `sudo apt remove
  mil-c` remove tudo.
- **Linux `.tar.gz`.** Para qualquer distribuição, sem administrador:

```
tar -xzf MIL-C-1.0.0-linux-x86_64.tar.gz
MIL-C-1.0.0/install.sh            # em ~/.local, com entrada no menu
MIL-C-1.0.0/install.sh --remove   # e para tirar
```

  O `install.sh` põe o programa em `~/.local/opt/MIL-C` e `milc` em `~/.local/bin`, e avisa se
  essa pasta não estiver no PATH. Sem instalar, `MIL-C-1.0.0/MIL-C/MIL-C` roda como está.
- As versões para Linux são feitas no Ubuntu 22.04 e rodam nele ou em qualquer sistema tão
  recente: o pacote usa a biblioteca C da máquina em que roda.

## A primeira abertura

Nenhum instalador é assinado, então a primeira abertura pede uma confirmação:

- **macOS**: clique com o botão direito (ou Control-clique) em **MIL-C** em Aplicativos, escolha
  **Abrir** e **Abrir** de novo. Se o macOS só oferecer mover para o lixo, abra **Ajustes do
  Sistema ▸ Privacidade e Segurança**, ache a linha sobre o MIL-C e escolha **Abrir Mesmo Assim**.
  Isso é pedido uma vez só.
- **Windows**: o SmartScreen diz que protegeu o computador. Escolha **Mais informações** e depois
  **Executar assim mesmo**.

O primeiro UMAP ou FlowSOM depois de instalar demora mais que os seguintes, porque o código
numérico é compilado uma vez e guardado na sua pasta de cache: medido no pacote do macOS, 16 s
para o primeiro UMAP e 5 s para cada um depois. [[troubleshooting]] trata de outros problemas na
primeira abertura.

## Conferindo o download

O `SHA256SUMS` traz o checksum de cada instalador, calculado sobre os arquivos que a página de
releases serve. Para conferir o seu:

```
shasum -a 256 -c SHA256SUMS --ignore-missing        # macOS
sha256sum -c SHA256SUMS --ignore-missing            # Linux
Get-FileHash MIL-C-1.0.0-windows-x64.msi            # Windows PowerShell: compare a olho
```

Todo instalador tem mais de 50 MB; um arquivo muito menor é um download quebrado.

## Conferindo o programa: selftest

Todo MIL-C instalado sabe se testar. O autoteste grava quatro arquivos FCS inventados numa pasta
temporária, abre-os como uma sessão, concatena, roda UMAP, tSNE e FlowSOM, grava uma figura em
SVG, PDF, PNG e TIFF, e depois faz o mesmo pela janela, sem mostrá-la. Ele imprime cada passo com
o tempo que levou e termina com `MIL-C 1.0.0: self-test passed`.

```
/Applications/MIL-C.app/Contents/MacOS/MIL-C selftest      # macOS
milc selftest                                              # Linux, instalado
MIL-C.exe selftest --log selftest.log                      # Windows, na pasta do programa
```

No Windows o programa não tem console, então `--log ARQUIVO` grava os passos, e qualquer falha,
num arquivo. `--no-window` pula a parte que precisa do Qt. Veja [[command-line]].

## Rodando a partir do código-fonte

O MIL-C é um pacote Python (`pyproject.toml`) e roda a partir de uma cópia do repositório em
Python 3.10 ou mais novo. Ele depende de numpy, scipy, scikit-learn, numba, umap-learn,
openTSNE, matplotlib, PyQt6 e pyqtgraph, que o pip instala:

```
git clone https://github.com/pedronobrez/mil-c
cd mil-c
pip install -e .
python -m milc                     # a janela
python -m milc gui caminho/da/pasta  # a janela, com uma pasta aberta
```

O `pip install -e .` também põe um comando `milc` no PATH, igual a `python -m milc`. No macOS,
um duplo clique em `MIL-C.command` no repositório também abre a janela.

Com o programa rodando, a próxima página é [[opening-data]].

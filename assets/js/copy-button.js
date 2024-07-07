document.addEventListener("DOMContentLoaded", function() {
    document.querySelectorAll('pre > code').forEach(function(codeBlock) {
        var button = document.createElement('button');
        button.className = 'copy-button';
        button.type = 'button';
        button.innerText = 'Copy';
        
        button.addEventListener('click', function() {
            var range = document.createRange();
            range.selectNode(codeBlock);
            window.getSelection().removeAllRanges(); // clear current selection
            window.getSelection().addRange(range); // to select text
            try {
                document.execCommand('copy');
                button.innerText = 'Copied';
                setTimeout(function() {
                    button.innerText = 'Copy';
                }, 2000);
            } catch (err) {
                console.log('Copy failed');
            }
            window.getSelection().removeAllRanges(); // to deselect
        });
        
        var pre = codeBlock.parentNode;
        if (pre.parentNode.classList.contains('highlight')) {
            pre.parentNode.insertBefore(button, pre);
        } else {
            pre.insertBefore(button, codeBlock);
        }
    });
});
